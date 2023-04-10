

#[derive(Clone,Debug)]
pub struct FastqRecord
{
    pub record: String,
    pub start_p7: usize,
    pub start_p5: usize,
    pub end_p7: usize,
    pub end_p5: usize,
    pub lane: Option<usize>
}

impl FastqRecord
{
    pub fn new ()
        -> Self
    {
        FastqRecord { record: String::new (), start_p7: 0, start_p5: 0, end_p7: 0, end_p5: 0, lane: None }
    }

    pub fn clear (&mut self)
    {
        self.record.clear ();
        self.start_p7 = 0;
        self.start_p5 = 0;
        self.end_p7 = 0;
        self.end_p5 = 0;
        self.lane = None;
    }
}

