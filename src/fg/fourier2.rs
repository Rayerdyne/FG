use super::spline::*;
use super::complex::*;
use std::f64::{self, consts::PI};
use super::fourier::CoeffsSet;

/// Integrates a cubic spline path and returns a set of 
/// Fourier coefficients. 
/// The path is described as 
/// (x(t), y(t)) = (sx(t), sy(t)), so that we integrate
/// sx(t) + j * sy(t) 
/// Description of method: 
/// https://www.overleaf.com/read/frhwfqrnkjjq 
#[allow(dead_code)]
pub fn compute_fourier_coeffs(sx: & Spline, sy: & Spline, n: usize) -> CoeffsSet {
    let (t_i, t_f) = (sx.start(), sx.end());
    assert_eq!(t_i, sy.start());
    assert_eq!(t_f, sy.end());
    assert_eq!(sx.num_parts(), sy.num_parts());

    // aside: T as a variable name makes rustc complain
    let period = t_f - t_i; 
    let omega_0_na = 2.0 * PI / period;
    let constants = Constants {
        omega_0_inv: FourTerms::new_pwr(1.0 / omega_0_na),
        changes: sx.changes()
    };

    let mut coeffs = CoeffsSet::new(n);

    coeffs.ppos[0] = Complex::zero();
    coeffs.nneg[0] = Complex::zero();
    for k in 1..n {
        coeffs.ppos[k] = compute_one( k as i32,   & sx, & sy, & constants);
        coeffs.nneg[k] = compute_one(-(k as i32), & sx, & sy, & constants);
    }

    coeffs
}

/// Computes the k-th fourier coefficient of sx(t) + j * sy(t). 
/// Achieves the sum over all the spline parts.
/// Output:  1/T * \hat f_k
fn compute_one(k: i32, sx: & Spline, sy: & Spline, constants:&  Constants) 
    -> Complex {
    
    let r = FourTerms::new( constants.omega_0_inv.na / k as f64,
        constants.omega_0_inv.sq / k.pow(2) as f64,
        constants.omega_0_inv.cu / k.pow(3) as f64,
        constants.omega_0_inv.fo / k.pow(4) as f64 );
    
    let mut x_k = Complex::zero();
    let mut y_k = Complex::zero();

    let mut t_i: ThreeTerms;
    let mut t_f: ThreeTerms = ThreeTerms::new_pwr(sx.start());
    for p in 0..sx.num_parts() {
        t_i = t_f;
        t_f = ThreeTerms::new_pwr(sx.changes()[p+1]);
        x_k += part_contribution(sx, p, & r, & t_i, & t_f);
        y_k += part_contribution(sy, p, & r, & t_i, & t_f);
    }

    (x_k + y_k.times_j()) / (2.0 * PI * constants.omega_0_inv.na)
}

/// Computes the contribution of the p-th part of the spline to the fourier
/// coefficient (of function s(t)) value, between t_i and t_f.
/// Output: \hat x_{k, p}
fn part_contribution(s: & Spline, p: usize, r: & FourTerms,
    t_i: & ThreeTerms, t_f: & ThreeTerms) -> Complex {

    let part = s.part(p);
    
    primitive(part, & t_f, & r) -  
    primitive(part, & t_i, & r)
}

/// Computes the value of the primitive (there is only one primitive...).
/// Output: F(k, a, b, c, d, t)
fn primitive(sp: SplinePart, t: & ThreeTerms, r: & FourTerms) -> Complex {
    let term_1 = r.na * (sp.a * t.cu + sp.b * t.sq + sp.c * t.na + sp.d);
    let term_2 = r.sq * (3.0 * sp.a * t.sq + 2.0 * sp.b * t.na + sp.c);
    let term_3 = r.cu * (6.0 * sp.a * t.na + 2.0 * sp.b);
    let term_4 = r.fo * (6.0 * sp.a);

    Complex {
        re: term_2 - term_4,
        im: term_1 - term_3
    } * Complex::expj(-t.na / r.na)
}

/// Holds 4 terms, which are powers of the na member 
#[derive(Debug)]
struct FourTerms {
    na: f64, // natural
    sq: f64, // squared
    cu: f64, // cubed
    fo: f64  // to the forth
}

impl FourTerms {
    fn new (na: f64, sq: f64, cu: f64, fo: f64) -> Self {
        Self {  na: na,     sq: sq,
                cu: cu,     fo: fo                  }
    }
    fn new_pwr (na: f64) -> Self {
        Self {  na: na,            sq: na.powi(2),
                cu: na.powi(3),    fo: na.powi(4)   }
    }
}

struct ThreeTerms {
    na: f64, // natural
    sq: f64, // squared
    cu: f64, // cubed
}

#[allow(dead_code)]
impl ThreeTerms {
    fn new (na: f64, sq: f64, cu: f64) -> Self {
        Self {  na: na,     sq: sq,
                cu: cu,          }
    }
    fn new_pwr (na: f64) -> Self {
        Self {  na: na,         sq: na.powi(2),
                cu: na.powi(3)   }
    }
}

/// Holds values that will remain constant to avoid unuseful recomputation.
#[derive(Debug)]
struct Constants {
    omega_0_inv: FourTerms,
    changes: Vec<f64>
}