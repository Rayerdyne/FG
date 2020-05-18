use super::spline::*;
use super::complex::*;
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

    // aside: T makes rustc complain
    let period = t_f - t_i; 
    let omega_0_na = 2.0 * PI / period;
    let constants = Constants {
        omega_0_inv: FourTerms::new_pwr(1.0 / omega_0_na),
        changes: sx.changes()
    };

    let mut coeffs = CoeffsSet::new(n);

    for k in 0..n {
        coeffs.ppos[k] = compute_one( k as i32,   & sx, & sy, & constants);
        coeffs.nneg[k] = compute_one(-(k as i32), & sx, & sy, & constants);
    }

    coeffs
}

/**
 * Computes the k-th fourier coefficient of sx(t) + j * sy(t). 
 * Achieves the sum over all the spline parts. 
 */
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

    x_k + y_k.times_j()
}

/**
 * Computes the contribution of the p-th part of the spline to the fourier
 * coefficient (of function s(t)) value.
 */
fn part_contribution(s: & Spline, p: usize, r: & FourTerms,
    t_i: & ThreeTerms, t_f: & ThreeTerms) -> Complex {

    let part = s.part(p);
    
    primitive(part, & t_i, & r) - primitive(part, & t_f, & r)
}

fn primitive(sp: SplinePart, t: & ThreeTerms, r: & FourTerms) -> Complex {
    let term_1 = r.na * (sp.a * t.cu + sp.b * t.sq + sp.c * t.na + sp.d);
    let term_2 = r.sq * (3.0 * sp.a * t.sq + 2.0 * sp.b * t.na + sp.c);
    let term_3 = r.cu * (6.0 * sp.a * t.na + 2.0 * sp.b);
    let term_4 = r.fo * (6.0 * sp.a);

    Complex {
        re: term_2 - term_4,
        im: term_1 - term_3
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

impl ThreeTerms {
    fn new (na: f64, sq: f64, cu: f64, fo: f64) -> Self {
        Self {  na: na,     sq: sq,
                cu: cu,          }
    }
    fn new_pwr (na: f64) -> Self {
        Self {  na: na,         sq: na.powi(2),
                cu: na.powi(3)   }
    }
}

/**
 * Holds values that will remain constant to avoid unuseful recomputation.
 */
#[derive(Debug)]
struct Constants {
    omega_0_inv: FourTerms,
    changes: Vec<f64>
}