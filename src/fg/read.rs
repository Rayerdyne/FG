use std::fs::File;
use std::io::Read;
use std::str::FromStr;
use std::fmt;
use std::num::ParseFloatError;

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

struct Point {
    x: f64,
    y: f64,
    t: f64,
}

#[allow(dead_code)]
pub struct PointsSet {
    xx: Vec<f64>,
    yy: Vec<f64>,
    tt: Vec<f64>,
}

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

impl FromStr for PointsSet {
    type Err = ReadingError;

    fn from_str(s: &str) -> Result<PointsSet, ReadingError> {
        let points_data: Vec<&str> = s.split('\n')
                                      .filter(|s| !s.is_empty())
                                      .collect();
        // println!("{:?}", points_data);
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

#[allow(dead_code)]
impl PointsSet {
    pub fn xx(&self)->Vec<f64> {    self.xx.clone()     }
    pub fn yy(&self)->Vec<f64> {    self.yy.clone()     }
    pub fn tt(&self)->Vec<f64> {    self.tt.clone()     }
}

impl std::convert::From<ParseFloatError> for ReadingError {
    fn from(e: ParseFloatError) -> ReadingError {
        ReadingError::ParseError(e)
    }
}

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


struct Complex {
    a: f64,
    b: f64,
}

impl FromStr for Complex {
    type Err = ReadingError;

    fn from_str(s: &str) -> Result<Complex, ReadingError> {
        let parts: Vec<&str> = s.trim_matches( |c| c == '(' || c == ')' )
                                .split(',').collect();
        let parsed_a: f64 = parts[0].trim().parse::<f64>()?;
        let parsed_b = parts[1].trim().parse::<f64>()?;
        Ok(Complex {
            a: parsed_a,
            b: parsed_b  })
    }
}

#[allow(dead_code)]
pub fn read_fourier_coeffs (filename: &str) 
        -> Result<(Vec<[f64; 2]>, Vec<[f64; 2]>), ReadingError> {
    let mut f = match File::open(filename) {
        Err(e) => return Err(ReadingError::FileStreamError(e, String::from(filename))),
        Ok(f) => f,
    };

    let mut data = String::new();
    if let  Err(e) = f.read_to_string(&mut data) {
        return Err(ReadingError::FileStreamError(e, String::from(filename)));
    }

    let lines: Vec<&str> =  data.split('\n')
                                .filter(|s| !s.is_empty())
                                .collect();

    let mut ppos = Vec::<[f64; 2]>::new();
    let mut nneg = Vec::<[f64; 2]>::new();
    for line in lines {
        let parts: Vec<&str> = line.split(' ').collect();
        if parts.len() < 2 {   return Err(ReadingError::IllFormedCoeffs)   }
        let cp = Complex::from_str(parts[0])?;
        let cn = Complex::from_str(parts[1])?;
        ppos.push([cp.a, cp.b]);
        nneg.push([cn.a, cn.b]);
        };

    Ok((ppos, nneg))
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