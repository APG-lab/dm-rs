

pub mod error
{
    use std::io;
    use thiserror::Error;

    #[derive(Error, Debug)]
    pub enum PublicError
    {
        #[error("IOError: {0}")]
        IOError (String)
    }

    impl From<io::Error> for PublicError
    {
        fn from (err: io::Error)
        -> PublicError
        {
            PublicError::IOError (err.to_string ())
        }
    }

}


