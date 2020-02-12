mod fourier;
mod fgif;
mod spline;
mod read;

extern crate clap;

use std::fmt;
use std::num::ParseIntError;
use clap::{Arg, App};

pub enum FgError {
    ReadingError(read::ReadingError),
    ParseArgError(ParseIntError),
    IoError(std::io::Error)
}

impl fmt::Display for FgError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FgError::ReadingError(e) => 
                { write!(f, "Read error: {}", e)           }
            FgError::ParseArgError(e) =>
                { write!(f, "Error parsing the provided arguments: {}", e)}
            FgError::IoError(e) =>
                { write!(f, "Error creating the output file: {}", e)}
        }
    }
}

impl std::convert::From<read::ReadingError> for FgError {
    fn from(e: read::ReadingError) -> FgError {
        FgError::ReadingError(e)
    }
}
impl std::convert::From<ParseIntError> for FgError {
    fn from(e: ParseIntError) -> FgError {
        FgError::ParseArgError(e)
    }
}
impl std::convert::From<std::io::Error> for FgError {
    fn from(e: std::io::Error) -> FgError {
        FgError::IoError(e)
    }
}

#[allow(dead_code)]
pub fn parse() -> Result<(), FgError> {
    
    let matches = App::new("fg")
            .version("0.1.0")
            .author("François Straet")
            .about("Drawings with Fourier series")
            .arg(Arg::with_name("input")
                .help("Sets the input file, containing the points of the drawing formatted as: `t: (x, y)`")
                .required(true)
                .index(1))
            .arg(Arg::with_name("output")
                .short("o")
                .long("output")
                .takes_value(true)
                .help("Sets the name of output file, `output.gif` if not provided"))
            .arg(Arg::with_name("fcolor")
                .short("f")
                .long("fcolor")
                .takes_value(true)
                .help("Sets the foreground color used in the output (hexcode)"))
            .arg(Arg::with_name("bcolor")
                .short("b")
                .long("bcolor")
                .takes_value(true)
                .help("Sets the background color used in the output (hexcode)"))
            .arg(Arg::with_name("width")
                .short("W")
                .long("gifwidth")
                .takes_value(true)
                .help("Sets the output's width"))
            .arg(Arg::with_name("height")
                .short("H")
                .long("gifheight")
                .takes_value(true)
                .help("Sets the output's height"))
            .arg(Arg::with_name("coeffs")
                .help("Ouputs drawing of custom Fourier coefficients in the input, which has to be formatted as `(c_k re, c_k im)&(c_-k re, c_-k im)`")
                .long("coeffs")
                .short("c"))
            .get_matches();

    let fcolor = matches.value_of("fcolor").unwrap_or("0x000000");
    let fc = color_from_hex(fcolor).unwrap();

    let bcolor = matches.value_of("bcolor").unwrap_or("0xFFFFFF");
    let bc = color_from_hex(bcolor).unwrap();

    let input = matches.value_of("input").unwrap();
    let output = matches.value_of("output").unwrap_or("output.gif");

    let sgw = matches.value_of("width").unwrap_or("300");
    let gw = sgw.parse::<usize>()?;
    let sgh = matches.value_of("height").unwrap_or("200");
    let gh = sgh.parse::<usize>()?;

    let coeffs_only = match matches.occurrences_of("coeffs"){
        0 => false,
        _ => true, };

    if coeffs_only {
        let coeffs = read::read_fourier_coeffs(input)?;
        println!("coeffs: \n{:?}", coeffs);

        fgif::draw_fourier_coeff(coeffs, output, gw, gh, 
            &[bc.0, bc.1, bc.2, fc.0, fc.1, fc.2])?;
    }
    else {
        let set = read::read_file(input)?;
        let ss = spline::interpolate_coords(vec![set.xx(), set.yy()], set.tt());
        let sx = ss[0].clone();
        let sy = ss[1].clone();
        println!("sx:\n {}", sx);
        println!("sy:\n {}", sy);    
    }
    Ok(())
}

#[allow(dead_code)]
pub fn test_gif(a: usize, b: usize, c: usize, d: usize) {
    fgif::gotest(a, b, c, d);
}

fn color_from_hex(s: &str) -> Result <(u8, u8, u8), ParseIntError> {
    let without_prefix = s.trim_start_matches("0x");
    println!("{}, then {}", s, without_prefix);
    let r = u8::from_str_radix(&without_prefix[0..2], 16)?;
    let g = u8::from_str_radix(&without_prefix[2..4], 16)?;
    let b = u8::from_str_radix(&without_prefix[4..6], 16)?;

    Ok((r, g, b))
}