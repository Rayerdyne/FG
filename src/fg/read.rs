use std::fs::File;
use std::io::Read;

use std::str::FromStr;
use std::num::ParseFloatError;

pub enum ReadingError {
    ParseError(ParseFloatError),
    FileStreamError(std::io::Error),
}

struct Point {
    x: f64,
    y: f64,
    t: f64,
}

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
    type Err = ParseFloatError;

    fn from_str(s: &str) -> Result<PointsSet, ParseFloatError> {
        let points_data: Vec<&str> = s.split('\n').collect();

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
                Err(e) => return Err(e)
            }
        }

        Ok(PointsSet{
            xx: parsed_xx,
            yy: parsed_yy, 
            tt: parsed_tt,
        })
    }
}

impl std::convert::From<std::io::Error> for ReadingError {
    fn from(e: std::io::Error) -> ReadingError {
        ReadingError::FileStreamError(e)
    }
}

impl std::convert::From<ParseFloatError> for ReadingError {
    fn from(e: ParseFloatError) -> ReadingError {
        ReadingError::ParseError(e)
    }
}

pub fn read_file (filename: &str) -> Result<PointsSet, ReadingError> {
    let mut f = File::open(filename)?;
    let mut data = String::new();
    f.read_to_string(&mut data)?;

    let set = PointsSet::from_str(&mut data)?;
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