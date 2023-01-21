# dm-rs
Demultiplex illumina sequence reads while allowing similar indices (i.e. edit distance &lt; 3)

This should dm a Novaseq run in under three hours. There is little optimisation so with a bit more work the bottlenecks could perhaps be reduced. Extending to directly read bcl files would be another useful thing. It's currently setup to use single threads for reading and writing and then multiple worker threads to process and distribute the records.

### Usage
This program reads Illumina reads from gzipped fastq files and requires the index sequences in the read name. It also requires a simple sample sheet tsv.

```bash
# cargo build --release
./target/release/dm \
	--sample-sheet-file-path <(echo -e "FOO\tACCTACTT\tGTCGCGGA\nBAR\tCCTGCCAA\tATAGTCAA") \
	--fastq-split-file-path examples/Undetermined_S0_L001_R1_001.test.fastq.gz \
	--output-name-pattern %sample%_R1.fastq.gz \
	--output-dir-path examples
# 0000	FOO	ACCTACTT	GTCGCGGA	1	0	1
# 0001	BAR	CCTGCCAA	ATAGTCAA	0	1	1


```


### Build
This program is written in [rust](https://www.rust-lang.org/) which uses the [cargo](https://doc.rust-lang.org/cargo/commands/cargo-build.html#compilation-options) build tool.

```bash
cargo build --release
```

