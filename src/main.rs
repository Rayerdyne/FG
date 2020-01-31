use std::env;
mod fg;
use std::num::ParseIntError;

fn main() -> Result<(), ParseIntError> {
    let args: Vec<String> = env::args().collect();

    if let Err(e) = fg::read(args) {
        println!("{}", e);
    };

    let reargs: Vec<String> = env::args().collect();
    let a = reargs[2].parse::<u32>()? as usize;
    let b = reargs[3].parse::<u32>()? as usize;
    let c = reargs[4].parse::<u32>()? as usize;
    let d = reargs[5].parse::<u32>()? as usize;

    fg::test_gif(a, b, c, d);
    Ok(())
}
