use std::fs::File;
use std::io::Read;
use std::str::FromStr;
use std::fmt;
use std::num::ParseFloatError;

use super::fourier::CoeffsSet;
use super::complex::Complex;

/// Represents error that could happen when reading files
pub enum ReadingError {
    ParseError(ParseFloatError),
    FileStreamError(std::io::Error, String),
    NotEnoughPoints,
    IllFormedCoeffs,
}

impl fmt::Display for ReadingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReadingError::ParseError(e) =>
                { write!(f, "Parsing error: {}", e)                        }
            ReadingError::FileStreamError(e, filename) => 
                { write!(f, "File stream error in `{}`: {}", filename, e)  }
            ReadingError::NotEnoughPoints =>
                { write!(f, "Not enough points !")                         }
            ReadingError::IllFormedCoeffs =>
                { write!(f, "Specified file is ill-formed !")              }
        }
    }
}

/// Holds the information needed for representing a point (a time and two
/// coordinates) plus the interpolation method that will be used BEFORE 
/// reaching that point, wether a line or a spline.
/// Members:
/// 
/// - `x`, `y`   : f64
/// - `t`       : f64
/// - `is_line` : bool
struct Point {
    x: f64,
    y: f64,
    t: f64,
    is_line: bool,
}

/// Holds a set of point to be interpolated
/// Each point has coordinates x, y and timestamp t. 
/// 
/// fromStr is implemented, for formatting:
/// (x, y, t)
/// (x, y, t)  ...
#[allow(dead_code)]
pub struct PointsSet {
    pub xx: Vec<f64>,
    pub yy: Vec<f64>,
    pub tt: Vec<f64>,
}

/// Parses a Point, the data should be formatted as:
/// `t: (x, y)`
/// Where `t` is the virtual time at which the point will be reached,
///       `x`, `y` are the coordinates of the point. Note that the axis will be
/// placed like in mathematics (origin at the center, `x` points to the right
/// and `y` to the left)
impl FromStr for Point {
    type Err = ParseFloatError;

    fn from_str(s: &str) -> Result<Point, ParseFloatError> {
        let parts: Vec<&str> = s.split(':')
                                .collect();
        let parsed_t = parts[0].trim()
                               .parse::<f64>()?;

        let coords: Vec<&str> = parts[1].trim()
                    .trim_matches( |c| c == '(' || c == ')' )
                    .split(',')
                    .collect();
        let parsed_x: f64 = coords[0].trim().parse::<f64>()?;
        let parsed_y: f64 = coords[1].trim().parse::<f64>()?;

        Ok(Point{x: parsed_x, y: parsed_y, t: parsed_t})
    } 
}

/// Parses a set of points, should be formatted with one point per line (cf 
/// Point::FromStr)
impl FromStr for PointsSet {
    type Err = ReadingError;

    fn from_str(s: &str) -> Result<PointsSet, ReadingError> {
        let points_data: Vec<&str> = s.split('\n')
                                      .filter(|s| !s.is_empty())
                                      .collect();
        if points_data.len() < 2 {
            return Err(ReadingError::NotEnoughPoints)
        }

        let mut parsed_xx = Vec::<f64>::new();
        let mut parsed_yy = Vec::<f64>::new();
        let mut parsed_tt = Vec::<f64>::new();

        for point_data in points_data {
            let p = Point::from_str(point_data);
            match p {
                Ok(point) => {
                    parsed_xx.push(point.x);
                    parsed_yy.push(point.y);
                    parsed_tt.push(point.t);
                }
                Err(e) => return Err(ReadingError::ParseError(e))
            }
        }

        Ok(PointsSet{
            xx: parsed_xx,
            yy: parsed_yy, 
            tt: parsed_tt,
        })
    }
}

impl std::convert::From<ParseFloatError> for ReadingError {
    fn from(e: ParseFloatError) -> ReadingError {
        ReadingError::ParseError(e)
    }
}

/// Parses a Complex formatted as
/// `(re, im)`
/// where `re` and `im` are doubles
impl FromStr for Complex {
    type Err = ReadingError;

    fn from_str(s: &str) -> Result<Complex, ReadingError> {
        let parts: Vec<&str> = s.trim_matches( |c| c == '(' || c == ')' )
                                .split(',').collect();
        let parsed_re: f64 = parts[0].trim().parse::<f64>()?;
        let parsed_im = parts[1].trim().parse::<f64>()?;
        Ok(Complex {
            re: parsed_re,
            im: parsed_im  })
    }
}

/// Parses a set of (Fourier) coefficient, i.e. an ensemble of complex numbers
/// The data should be n lines formatted as:
/// `<Complex_x> & <Complex_y>`
/// Where the line number indicates the rank of the coefficient (i-th line -> 
/// i-th harmonic), `Complex_x` is the coefficient for the x(t) function and 
/// `Complex_y` for the y(t) function.
impl FromStr for CoeffsSet {
    type Err = ReadingError;

    fn from_str(s: &str) -> Result<CoeffsSet, ReadingError> {
        let lines: Vec<&str> =  s.split('\n')
        .filter(|s| !s.is_empty())
        .collect();

        let mut ppos = Vec::<Complex>::new();
        let mut nneg = Vec::<Complex>::new();
        for line in lines {
            let parts: Vec<&str> =  line.split('&')
                        .filter(|s| !s.is_empty())
                        .collect();
            if parts.len() < 2 {   return Err(ReadingError::IllFormedCoeffs)   }
            let cp = Complex::from_str(parts[0])?;
            let cn = Complex::from_str(parts[1])?;
            ppos.push(cp);
            nneg.push(cn);    
        };
        Ok(CoeffsSet{ppos: ppos, nneg: nneg})
    }
}

/// Reads a file and return the set of points it contains. 
///  RETURN       Result<PointsSet, ReadingError>
#[allow(dead_code)]
pub fn read_file(filename: & str) -> Result<PointsSet, ReadingError> {
    let mut f = match File::open(filename) {
        Err(e) => return Err(ReadingError::FileStreamError(e, String::from(filename))),
        Ok(f) => f,
    };

    let mut data = String::new();
    if let  Err(e) = f.read_to_string(&mut data) {
        return Err(ReadingError::FileStreamError(e, String::from(filename)));
    }

    let set = PointsSet::from_str(&mut data)?;
    Ok(set)
}

/// Reads Fourier coefficients written in the file named filename. 
/// RETURN           Result<CoeffsSet, ReadingError>
#[allow(dead_code)]
pub fn read_fourier_coeffs (filename: &str) -> Result<CoeffsSet, ReadingError> {
    let mut f = match File::open(filename) {
        Err(e) => return Err(ReadingError::FileStreamError(e, String::from(filename))),
        Ok(f) => f,
    };

    let mut data = String::new();
    if let  Err(e) = f.read_to_string(&mut data) {
        return Err(ReadingError::FileStreamError(e, String::from(filename)));
    }

    let set = CoeffsSet::from_str(&mut data)?;
    Ok(set)
}


/*
Parsing floats accepts following formattings:
    '3.14'
    '-3.14'
    '2.5E10', or equivalently, '2.5e10'
    '2.5E-10'
    '5.'
    '.5', or, equivalently, '0.5'
    'inf', '-inf', 'NaN'

    Copied-pasted from https://doc.rust-lang.org/std/str/trait.FromStr.html
*/