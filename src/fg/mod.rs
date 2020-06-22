mod complex;
mod fourier;
mod fgif;
mod spline;
mod read;

extern crate clap;

use std::fmt;
use std::num::ParseIntError;
use clap::{Arg, App};
use std::f64::consts::PI;

/// Error type returned by `parse()` function,
/// represents what's could go wrong.
/// 
/// fmt::Display is logically implemented =)
#[allow(dead_code)]
pub enum FgError {
    ReadingError(read::ReadingError),
    IoError(std::io::Error)
}

impl fmt::Display for FgError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FgError::ReadingError(e) => 
                { write!(f, "Read error: {}", e)           }
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
impl std::convert::From<std::io::Error> for FgError {
    fn from(e: std::io::Error) -> FgError {
        FgError::IoError(e)
    }
}

const STD: u8 = 1;
const COEFFS_ONLY: u8 = 2;
const SPLINE: u8 = 3;

const DEF_HEIGHT: usize = 200;
const DEF_WIDTH: usize = 300;
const DEF_N_STEPS: usize = 200;
const DEF_N_COEFFS: usize = 5;
/// Parses arguments provided to the program, and process to execution
#[allow(dead_code)]
pub fn parse() -> Result<(), FgError> {
    
    let matches = app_args();

    let fc = get_color(& matches, "fcolor", "0x000000");
    let bc = get_color(& matches, "bcolor", "0xFFFFFF");

    let input = matches.value_of("input").unwrap();
    let output = matches.value_of("output").unwrap_or("output.gif");

    let gw = get_value(& matches, "width", DEF_WIDTH);
    let gh = get_value(& matches, "width", DEF_HEIGHT);
    let n_steps = get_value(& matches, "width", DEF_N_STEPS);
    let n_coeffs = get_value(& matches, "n_coeffs", DEF_N_COEFFS) + 1;

    let ctype = match matches.value_of("test").unwrap_or("std") {
        "std" => STD,
        "coeffs" => COEFFS_ONLY,
        "spline" => SPLINE,
        _ => STD,
    };

    if ctype == COEFFS_ONLY {
        let coeffs = read::read_fourier_coeffs(input)?;
        println!("coeffs: \n{}", coeffs);

        fgif::draw_fourier_coeff(coeffs, output, gw, gh, (0.0, 2.0*PI), n_steps,
            &[bc.0, bc.1, bc.2, fc.0, fc.1, fc.2])?;
    }
    else if ctype == SPLINE {
        let set = read::read_file(input)?;
        let ss = spline::interpolate_coords(vec![set.xx, set.yy], set.tt);
        let sx = ss[0].clone();
        let sy = ss[1].clone();

        fgif::draw_spline(sx, sy, output, gw, gh, n_steps, 
            &[bc.0, bc.1, bc.2, fc.0, fc.1, fc.2])?;
    }
    else {
        let set = read::read_file(input)?;
        let ss = spline::interpolate_coords(vec![set.xx, set.yy], set.tt);
        let sx = ss[0].clone();
        let sy = ss[1].clone();
        
        let coeffs = fourier::compute_fourier_coeffs(& sx, & sy, n_coeffs);

        println!("{}", coeffs);
        fgif::draw_fourier_coeff(coeffs, output, gw, gh, (sx.start(), sx.end()),
            n_steps, &[bc.0, bc.1, bc.2, fc.0, fc.1, fc.2])?;  
    }
    println!("Wrote {} frames in {} ({}, {}), with {} coeffs", n_steps, output,
        gw, gh, n_coeffs);
    Ok(())
}

fn app_args() -> clap::ArgMatches<'static> {
    App::new("fg")
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
        .arg(Arg::with_name("n_steps")
                .short("n")
                .long("n-steps")
                .takes_value(true)
                .help("Sets the numbers of frames of the output."))
        .arg(Arg::with_name("n_coeffs")
            .short("c")
            .long("n-coeffs")
            .takes_value(true)
            .help("Sets Fourier coefficients computed and used."))
        // .arg(Arg::with_name("coeffs")
        //     .help("Ouputs drawing of custom Fourier coefficients in the input, which has to be formatted as\n \
        //             `(Re(c_k),Im(c_k))&(Re(c_-k) , Im(c_-k))`")
        //     .long("coeffs")
        //     .short("c"))
        .arg(Arg::with_name("test")
            .takes_value(true)
            .possible_values(&["coeffs", "spline"])
            .help("Coeffs: uputs drawing of custom Fourier coefficients in the input, which has to be formatted as\n \
                    `(Re(c_k),Im(c_k))&(Re(c_-k) , Im(c_-k))` \n
                    Spline: only draws the spline.")
            .long("type")
            .short("t"))
        .get_matches()
}

#[allow(dead_code)]
pub fn test_gif(a: usize, b: usize, c: usize, d: usize) {
    fgif::gotest(a, b, c, d);
}

/// Get a color desribed in argument, default value if not present.
fn get_color(matches: & clap::ArgMatches, name: & str, def: & str) 
    -> (u8, u8, u8) {
    let s = matches.value_of(name).unwrap_or(def);
    color_from_hex(s).unwrap()
}

/* Parses a color written in hexcode (0xrrggbb) to a tuple (r, g, b) */
fn color_from_hex(s: &str) -> Result <(u8, u8, u8), ParseIntError> {
    let without_prefix = s.trim_start_matches("0x");
    // println!("{}, then {}", s, without_prefix);
    let r = u8::from_str_radix(&without_prefix[0..2], 16)?;
    let g = u8::from_str_radix(&without_prefix[2..4], 16)?;
    let b = u8::from_str_radix(&without_prefix[4..6], 16)?;

    Ok((r, g, b))
}

fn get_value(matches: & clap::ArgMatches, name: & str, def: usize) -> usize {
    let def_str = def.to_string();
    let s = matches.value_of(name).unwrap_or(& def_str);
    s.parse::<usize>().unwrap_or(def)
}