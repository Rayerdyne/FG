use super::spline::*;
use std::f64::{self, consts::PI};
use std::fmt;

/** Integrates a cubic spline path and returns a set of 
 * Fourier coefficients. 
 * The path is described as 
 * (x(t), y(t)) = (sx(t), sy(t)), so that we integrate
 * sx(t) + j * sy(t) 
 * Description of method: 
 * https://www.overleaf.com/read/frhwfqrnkjjq 
 * */
#[allow(dead_code)]
pub fn compute_fourier_coeffs(sx: Spline, sy: Spline, n: usize) -> CoeffsSet {
    let (t_i, t_f) = (sx.start(), sx.end());
    assert_eq!(t_i, sy.start());
    assert_eq!(t_f, sy.end());
    assert_eq!(sx.num_parts(), sy.num_parts());

    let T = t_f - t_i;
    let omega_0_na = 2.0 * PI / T;
    let omega_0_inv = QuadTerm::new_pwr(1 / omega_0_na);

    let mut coeffs = CoeffsSet::new(n);

    for k in 0..n {
        coeffs.ppos[k] = compute_one( k, sx, sy, omega_0_inv);
        coeffs.nneg[k] = compute_one(-k, sx, sy, omega_0_inv);
    }

    coeffs
}

/**
 * Computes the k-th fourier coefficient of sx(t) + j * sy(t). 
 * Achieves the sum over all the spline parts. 
 */
fn compute_one(k: usize, sx: Spline, sy: Spline, omega_0_inv: QuadTerm) -> Complex {


}

/** Holds a complex number */
#[derive(Debug, Clone, Copy)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl std::ops::AddAssign for Complex {
    fn add_assign(&mut self, c: Complex) {  self.re += c.re;
                                            self.im += c.im;    }
}
impl std::ops::DivAssign<f64> for Complex {
    fn div_assign(&mut self, div: f64) {    self.re /= div;
                                            self.im /= div;     }
}

impl Complex {
    fn times_i(&mut self) -> Complex {
        let temp = self.re;
        self.re = -self.im;
        self.im = temp;
        *self
    }
    fn zero() -> Complex {
        Complex {   re: 0.0,
                    im: 0.0     }
    }
}

/** Holds a set of Fourier coefficients. */
#[derive(Debug)]
pub struct CoeffsSet {
    pub ppos: Vec<Complex>,
    pub nneg: Vec<Complex>,
    //doubled character because vector
}

impl CoeffsSet {
    fn new(n: usize) -> CoeffsSet {
        CoeffsSet {     ppos: vec![Complex::zero(); n],
                        nneg: vec![Complex::zero(); n],     }
    }
}

impl fmt::Display for CoeffsSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        for i in 0..self.ppos.len() {
            s.push_str(&format!("{}:\t+{:?}\n\t-{:?}\n", i, self.ppos[i], self.nneg[i]));
        }
        write!(f, "{}", s)
    }
}

/** Holds 4 terms, which are powers of the na member */
#[derive(Debug)]
struct QuadTerm {
    na: f64, // natural
    sq: f64, // squared
    cu: f64, // cubed
    fo: f64  // to the forth
}

impl QuadTerm {
    fn new (na: f64, sq: f64, cu: f64, fo: f64) {
        QuadTerm {  na: na,     sq: sq,
                    cu: cu,     fo: fo      }
    }
    fn new_pwr (na: f64) {
        QuadTerm {  na: na,         sq: na.powi(2),
                    cu: na.powi(3), fo: na.powi(4)    }
    }
}

#[derive(Debug)]
struct Constants {
    omega_0_inv: QuadTerm,

}