use std::env;
mod fg;

fn main(){
    let args: Vec<String> = env::args().collect();

    if let Err(e) = fg::read(args) {
        println!("{}", e);
    };
}
