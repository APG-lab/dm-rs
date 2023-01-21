
use clap::Parser;
use crossbeam_channel;
use flate2::Compression;
use flate2::write::GzEncoder;
use log::{debug,error};
use std::collections;
use std::io::Write;
use std::iter;
use std::fs;
use std::process;
use std::time;
use std::sync;
use tokio;

mod dm;
mod file;
mod helper;

#[derive(Debug,Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli
{
    #[arg(long, help="Split fastq files to dm",num_args=1..)]
    fastq_split_file_path: Vec<String>,

    #[arg(long, help = "Sample sheet file path")]
    sample_sheet_file_path: String,

    #[arg(long, help = "Number of mismatches allowed", default_value="1")]
    n_mismatches: usize,

    #[arg(long, help = "Number of worker threads", default_value="3")]
    n_threads: usize,

    #[arg(long, help = "Path to output directory")]
    output_dir_path: String,

    #[arg(long, help = "Output name pattern. %sample% will be substituted. Should end in '.fastq.gz'", default_value="%sample%.fastq.gz")]
    output_name_pattern: String
}

fn main () {
    env_logger::init ();

    let reads_in = sync::Arc::new (sync::atomic::AtomicUsize::new (0));
    let reads_in_progress = reads_in.clone ();

    let dbg_reads_proc = sync::Arc::new (sync::atomic::AtomicUsize::new (0));
    let dbg_reads_proc_end = dbg_reads_proc.clone ();

    let progress_ms = 10_000;
    let chan_bound = 100_000;

    debug! ("hello");
    let args = Cli::parse ();

    //debug! ("args: {:?}",args);
    let sample_data = file::read_tsv (args.sample_sheet_file_path).expect ("Failed to read sample sheet");
    //debug! ("sample_data: {:?}",sample_data);

    if args.fastq_split_file_path.len () == 0
    {
        eprintln! ("Please supply at least one fastq file to demultiplex");
        process::exit (1);
    }
    else
    {
        let rt = tokio::runtime::Builder::new_multi_thread ()
        .worker_threads (args.n_threads)
        .enable_all ()
        .build ()
        .unwrap ();

        let (shutdown_send, mut shutdown_recv) = tokio::sync::mpsc::unbounded_channel::<()> ();
        let (tx, _) = tokio::sync::broadcast::channel::<()> (args.n_threads + 4);
        let tx_cancel = tx.clone ();
        /*
         * rx_stats was going out of scope before the stats were sent. I could use
         * crossbeam_channel and clone the receiver to keep it alive, but better to move
         * this to a blocking thread before the program exits?
         */
        let (tx_stats, mut rx_stats) = tokio::sync::mpsc::unbounded_channel::<((String, String), usize, usize)> ();
        let (tx_worker, rx_worker) = crossbeam_channel::bounded::<dm::FastqRecord> (args.n_threads);

        // create dm channels
        let (mut tx_sample_data_dm, mut rx_sample_data_dm) = sample_data.iter ().enumerate ().fold ( ( Vec::new (), Vec:: new () ), |mut acc, (i, item)| {
            let (tx_sample, rx_sample) = tokio::sync::mpsc::channel::<dm::FastqRecord> (chan_bound);
            let output_name_resolved = args.output_name_pattern.replace ("%sample%", &item[0]);
            let output_file_path = format! ("{}/{}", args.output_dir_path, output_name_resolved);
            let output_file = GzEncoder::new (fs::File::create (output_file_path.clone ()).expect (&format! ("Unable to create file {}", output_file_path.clone ())), Compression::default ());
            acc.0.push ( ( item, tx_sample ) );
            acc.1.push ( ( rx_sample, output_file, i, output_file_path) );
            acc
        });
        // Ensure our index pairs are sorted so we can use binary search
        tx_sample_data_dm.sort_by (|a, b| (&a.0[1], &a.0[2]).partial_cmp ( &( &b.0[1], &b.0[2] ) ).unwrap());

        for (i, (item, _) ) in tx_sample_data_dm.iter ().enumerate ()
        {
            debug! ("{:04} {}+{}", i,  item[1], item[2]);
        }

        // spawn workers
        for i in 1..args.n_threads
        {
            let mut index_lookup = tx_sample_data_dm.iter ().fold ( ( Vec::new (), Vec::new (), Vec::new (), Vec::new () ), |mut acc, item| {
                acc.0.push ( ( item.0[1].clone (), item.0[2].clone () ) );
                acc.1.push ( item.1.clone () );
                acc.2.push ( 0 );
                acc.3.push ( 0 );
                acc
            });
            //debug! ("index_lookup: {:?}", index_lookup);
            let worker_tx_stats = tx_stats.clone ();
            let worker_rx_worker = rx_worker.clone ();
            let mut worker_rx_cancel = tx.subscribe ();
            let worker_dbg_reads_proc = dbg_reads_proc.clone ();
            rt.spawn_blocking (move || {
                loop
                {
                    match worker_rx_cancel.try_recv ()
                    {
                        Ok (_) => {
                            debug! ("{} abort worker", i);
                            break;
                        },
                        Err (tokio::sync::broadcast::error::TryRecvError::Closed) => {
                            debug! ("{} worker closed", i);
                            break;
                        },
                        _ => {}
                    }
                    match worker_rx_worker.recv ()
                    {
                        Ok (msg) => {
                            //debug! ("{} dm got msg: {:?}", i, msg);
                            worker_dbg_reads_proc.fetch_add (1,sync::atomic::Ordering::Relaxed);
                            if let ( Some (inda), Some (indb) ) = ( &msg.record.get (msg.start_p7..msg.end_p7), &msg.record.get (msg.start_p5..msg.end_p5) )
                            {
                                match index_lookup.0.binary_search ( &(inda.to_string (), indb.to_string ()) )
                                {
                                    Ok (sample_index) => {
                                        match index_lookup.1[sample_index].blocking_send (msg)
                                        {
                                            Ok (_) => {
                                                index_lookup.2[sample_index] += 1;
                                            },
                                            Err (e) => {
                                                debug! ("Failed to send to sample receiver: {:?}", e);
                                            }
                                        }
                                    },
                                    Err (_) => {
                                        if args.n_mismatches > 0
                                        {
                                            // attempt to find with 1 mismatch
                                            for (sample_index, (sinda, sindb) ) in index_lookup.0.iter ().enumerate ()
                                            {
                                                if iter::zip (sinda.chars (), inda.chars ()).filter (|(x,y)| x != y).count () <= args.n_mismatches && iter::zip (sindb.chars (), indb.chars ()).filter (|(x,y)| x != y).count () <= args.n_mismatches
                                                {
                                                    match index_lookup.1[sample_index].blocking_send (msg)
                                                    {
                                                        Ok (_) => {
                                                            index_lookup.3[sample_index] += 1;
                                                        },
                                                        Err (e) => {
                                                            debug! ("Failed to send to sample receiver: {:?}", e);
                                                        }
                                                    }
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                        Err (_e) => {
                            //debug! ("{} recv err: {:?}", i, e);
                            break;
                        }
                    }
                }
                drop (worker_rx_worker);
                for sndr in index_lookup.1.iter_mut ()
                {
                    drop (sndr);
                }
                debug! ("worker {} shutdown", i);
                for ( (inda, indb), (n_perfect, n_mm) ) in iter::zip (index_lookup.0, iter::zip (index_lookup.2, index_lookup.3))
                {
                    //println! ("{}\t{}\t{}\t{}\t{}\t{}", sample_name, inda, indb, n_perfect, n_mm, n_perfect + n_mm);
                    match worker_tx_stats.send ( ( (inda, indb), n_perfect, n_mm) )
                    {
                        Ok (_) => {},
                        Err (e) => {
                            debug! ("{} Failed to send stats {:?}", i, e);
                        }
                    }
                }
            });
        }
        // drop the original worker_receiver, or our fastq read loop won't cancel
        drop (rx_worker);
        // We need to drop the original senders
        for (_, sndr) in tx_sample_data_dm.iter_mut ()
        {
            drop (sndr);
        }
        debug! ("dropped original tx_sample senders");

        rt.spawn (async move {
            let mut progress_interval = tokio::time::interval (time::Duration::from_millis (progress_ms));
            let mut last_reads_in = 0;
            loop {
                // returns Result of the evaluated handler
                match tokio::select! {
                    _ = tokio::signal::ctrl_c () => { tx_cancel.send (()).ok ();Err ("Got crl-c") },
                    _ = shutdown_recv.recv () => { tx_cancel.send (()).ok ();Err ("All senders out of scope") },
                    _ = progress_interval.tick () => {
                        let current_reads_in = reads_in_progress.load (sync::atomic::Ordering::Relaxed);
                        let rate = ((current_reads_in - last_reads_in) as f64 * (1000 as f64/progress_ms as f64)) as usize;
                        debug! ("in: {:012} rate: {:09}/s", current_reads_in, rate);
                        last_reads_in = current_reads_in;
                        Ok (())
                    },
                }
                {
                    Ok (_) => (),
                    Err (err) => { debug! ("select loop shutdown {}", err);break; }
                }
            }
        });

        // Listen to sample channels and write to files
        rt.spawn (async move {
            let mut closed = (0..rx_sample_data_dm.len ()).collect::<collections::HashSet<usize>> ();
            loop
            {
                for (chan, wtr, i, fp) in rx_sample_data_dm.iter_mut ()
                {
                    for _ in 0..chan_bound
                    {
                        match chan.try_recv ()
                        {
                            Err (tokio::sync::mpsc::error::TryRecvError::Empty) => {
                                break;
                            },
                            Ok (msg) => {
                                match wtr.write_all (msg.record.as_bytes ())
                                {
                                    Ok (_) => {},
                                    Err (e) => {
                                        error! ("Failed to write to '{}' reason: {:?}", fp, e);
                                    }
                                }
                            },
                            Err (tokio::sync::mpsc::error::TryRecvError::Disconnected) => {
                                // No more records will be sent
                                match wtr.try_finish ()
                                {
                                    Ok (_) => {},
                                    Err (e) => {
                                        error! ("Failed final write to '{}' reason: {:?}", fp, e);
                                    }
                                }
                                closed.remove (i);
                            }
                        }
                    }
                }
                if closed.len () == 0
                {
                    debug! ("shutdown write loop");
                    break;
                }
            }
        });

        //debug! ("args.fastq_split_file_path: {:?}", args.fastq_split_file_path);
        //file::read_fastq (args.fastq_split_file_path[0].clone (), dm::process_fastq_line);

        let process_fastq_line = |line: &str, line_number: &usize, records_read: &mut usize, fastq_record: &mut dm::FastqRecord| -> bool {
            if line_number % 4 == 0 && *line_number > 0
            {
                if line_number % 100_000 == 0
                {
                    reads_in.fetch_add (25_000, sync::atomic::Ordering::Relaxed);
                }
                // find indices
                if let Some (eofl) = fastq_record.record.find ("\n")
                {
                    //debug! ("found eofl: {}", eofl);
                    if let Some (solocs) = fastq_record.record[0..eofl].rfind (":")
                    {
                        if let Some (sind) = fastq_record.record[solocs..eofl].find ("+")
                        {
                            if let (Some (_inda), Some (_indb)) = ( fastq_record.record.get ((solocs+1)..(solocs+sind)), fastq_record.record.get ((solocs+1+sind)..eofl) )
                            {
                                //debug! ("inda: {}", &fastq_record.record[solocs..(solocs+sind)]);
                                //debug! ("inda: {} indb: {}", inda, indb);
                                fastq_record.start_p7 = solocs+1;
                                fastq_record.start_p5 = solocs+1+sind;
                                fastq_record.end_p7 = solocs+sind;
                                fastq_record.end_p5 = eofl;
                            }
                            else
                            {
                                panic! ("Failed to obtain indices from calculated offsets in fastq record header");
                            }
                        }
                        else
                        {
                            panic! ("Failed to find index sep '+' in fastq record header");
                        }
                    }
                    else
                    {
                        panic! ("Failed to find indices sep ':' in fastq record header");
                    }
                }
                else
                {
                    panic! ("Failed to find first line of fastq record");
                }
                *records_read += 1;
                match tx_worker.send (fastq_record.clone ())
                {
                    Ok (_) => {},
                    Err (e) => {
                        debug! ("Failed to send fastq_record: {:?}", e);
                        return false;
                    }
                }
                fastq_record.clear ();
            }
            fastq_record.record.push_str (line);
            fastq_record.record.push_str ("\n");
            return true;
        };
        let mut records_read = 0;
        for fastq_file_path in args.fastq_split_file_path
        {
            if !file::read_fastq (&mut records_read, fastq_file_path, process_fastq_line).expect ("File was not read properly")
            {
                debug! ("aborting fastq file reading");
            }
        }
        debug! ("total records read: {}", records_read);

        drop (shutdown_send);
        debug! ("gather stats");
        rt.spawn_blocking (move || {
            let mut stats = collections::HashMap::<(String,String), (usize,usize)>::new (); 
            while let Some ( (inds, n_perfect, n_mm ) ) = rx_stats.blocking_recv ()
            {
                let stat = stats.entry (inds).or_insert ( ( 0, 0) );
                stat.0 += n_perfect;
                stat.1 += n_mm;
            }
            //debug! ("stats: {:?}", stats);
            for (i, item) in sample_data.iter ().enumerate ()
            {
                if let Some (stat) = stats.get ( &(item[1].clone (), item[2].clone ()) )
                {
                    println! ("{:04}\t{}\t{}\t{}\t{}\t{}\t{}", i, item[0], item[1], item[2], stat.0, stat.1, stat.0 + stat.1);
                }
                else
                {
                    println! ("{:04}\t{}\t{}\t{}\t?\t?\t?", i, item[0], item[1], item[2]);
                }
            }
        });
        debug! ("reads_proc: {}", dbg_reads_proc_end.load (sync::atomic::Ordering::Relaxed));
        debug! ("Server shutting down");
    }
}


