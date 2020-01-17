mod fourier;
mod gif;
mod spline;

pub fn tell() {
    println!("In tell function in /fg/mod.rs" );
    fourier::fourier();
    gif::gif();
    spline::spline();
}