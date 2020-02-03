mod fg;

fn main() {

    if let Err(e) = fg::parse() {
        println!("{}", e);
    };

}
