mod fourier;
mod gif;
mod spline;
mod read;

use std::fmt;

pub enum FgError {
    ArgError,
    ReadingError(read::ReadingError),
}

impl fmt::Display for FgError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FgError::ArgError =>
                { write!(f, "Arguments are not valid !")   }
            FgError::ReadingError(e) => 
                { write!(f, "Read error: {}", e)           }
        }
    }
}

impl std::convert::From<read::ReadingError> for FgError {
    fn from(e: read::ReadingError) -> FgError {
        FgError::ReadingError(e)
    }
}

#[allow(dead_code)]
pub fn read(args: Vec<String>) -> Result<(), FgError> {
    if args.len() <= 1 {    return Err(FgError::ArgError)   }
    
    let set = read::read_file(&args[1])?;
    let sx = spline::interpolate(set.get_xx(), set.get_tt());
    let sy = spline::interpolate(set.get_yy(), set.get_tt());
    println!("sx: {}", sx);
    println!("sy: {}", sy);
    Ok(())
}