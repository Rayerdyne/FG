mod fourier;
mod gif;
mod spline;
mod read;

pub fn tell() {
    println!("In tell function in /fg/mod.rs" );
    fourier::hello_fourier();
    gif::hello_gif();
    spline::hello_spline();
}
