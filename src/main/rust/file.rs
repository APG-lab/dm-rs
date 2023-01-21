
use crate::dm;
use crate::helper;
use flate2::read::GzDecoder;
//use log::debug;
use std::io;
use std::fs;
use std::io::prelude::*;

pub fn read_tsv (tsv_file_path: String)
    -> Result<Vec<Vec<String>>, helper::error::PublicError>
{
    let mut f = fs::File::open (tsv_file_path)?;

    let mut contents = String::new ();
    f.read_to_string (&mut contents)?;

    //println! ("read {:?}", contents);

    let result: Vec<Vec<String>> = contents.trim_end_matches ("\n").split ("\n").map (| line | -> Vec<String> {
        let ldata: Vec<String> = line.split ("\t").map (String::from).collect ();
        ldata
    }).collect ();

    //println! ("result {:?}", result);
    Ok (result)
}


pub fn read_fastq<F> (records_read: &mut usize, fastq_file_path: String, mut process_line: F)
    -> Result<bool, helper::error::PublicError>
    where
        F: FnMut (&str, &usize, &mut usize, &mut dm::FastqRecord) -> bool
{

    let mut reader = io::BufReader::new (GzDecoder::new (fs::File::open (fastq_file_path)?));
    let mut buffer = String::new ();
    let mut line_count: usize = 0;
    let mut fastq_record = dm::FastqRecord::new ();
    let mut read_ok = true;
    while let Some(n) = reader.read_line(&mut buffer).ok ()
    {
        //debug! ("n: {}",n);
        //debug! ("{}", buffer.trim ());
        if n == 0
        {
            break;
        }
        if process_line (buffer.trim (), &line_count, records_read, &mut fastq_record)
        {
            buffer.clear ();
            line_count += 1;
        }
        else
        {
            read_ok = false;
            break;
        }
    }
    if read_ok && process_line (buffer.trim (), &line_count, records_read, &mut fastq_record)
    {
        Ok (true)
    }
    else
    {
        Ok (false)
    }
}
