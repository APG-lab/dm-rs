# dm-rs
Demultiplex illumina sequence reads while allowing similar indices (i.e. edit distance &lt; 3)

This should dm a Novaseq run in under three hours. There is little optimisation so with a bit more work the bottlenecks could perhaps be reduced. Extending to directly read bcl files would be another useful thing. It's currently setup to use single threads for reading and writing and then multiple worker threads to process and distribute the records.
