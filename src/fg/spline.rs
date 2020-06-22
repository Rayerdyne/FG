extern crate nalgebra as na;

use na::{DMatrix, DVector};
use std::fmt;

/// Represents a cubic spline.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct Spline {
    /// Contains each cubic polynimial coefficients, i.e. a, b, c and d.
    parts: Vec<SplinePart>,
    /// Contains times at which the polynomial that is equal to the spline
    /// changes.
    changes: Vec<f64>,
    /// For the iterator...
    current: usize,
    /// Time at which the spline start, i.e. `spline(t < start) = 0`
    start: f64,
    /// Time at which the spline ends, i.e. `spline(t > end) = 0`
    end: f64,
}

/// Represents a spline part, i.e. a cubic polynomial. 
/// Thus, it holds its coefficients a, b, c, d
#[allow(dead_code)]
#[derive(Copy, Clone, Debug)]
pub struct SplinePart {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
}

/// Returns the matrix that will have to be solved when
/// computing a spline interpolating points whose timestamps
/// are in tt vector 
fn matrix_for (tt: &Vec<f64>, mm: &Vec<bool>) -> DMatrix<f64> {
    let n = tt.len();
    let mut a: DMatrix<f64> = DMatrix::zeros(4*(n-1), 4*(n-1));

    let mut t1_squared: f64 = tt[0] * tt[0];
    let mut t1_cubed:   f64 = t1_squared * tt[0];
    let mut t2_squared: f64 = tt[1] * tt[1];
    let mut t2_cubed:   f64 = t2_squared * tt[1];

    for i in 0..(n-1) {
        if i >= 1 {
            t1_squared = t2_squared;
            t1_cubed = t2_cubed;
        }
        t2_squared = tt[i+1] * tt[i+1];
        t2_cubed = t2_squared * tt[i+1];

        a[(4*i, 4*i)]   = t1_cubed;
        a[(4*i, 4*i+1)] = t1_squared;
        a[(4*i, 4*i+2)] = tt[i];
        a[(4*i, 4*i+3)] = 1.0;

        a[(4*i+1, 4*i)]   = t2_cubed;
        a[(4*i+1, 4*i+1)] = t2_squared;
        a[(4*i+1, 4*i+2)] = tt[i+1];
        a[(4*i+1, 4*i+3)] = 1.0;
    
        if i < n-2 {
            add_cont_equations(&mut a, i, tt[i+1], mm[i]);
        }
    }

    /* Starting & ending slopes are 0 : two missing eqations */
    if mm[0] {
        a[(4*n-6, 0)] = 1.0;
    } else {
        a[(4*n-6, 0)] = tt[0] * tt[0];
        a[(4*n-6, 1)] = tt[0];
        a[(4*n-6, 2)] = 1.0;
    }
    
    if mm[n-1] {
        a[(4*n-5, 4*n-8)] = 1.0;
    } else {
        a[(4*n-5, 4*n-8)] = tt[n-1] * tt[n-1];
        a[(4*n-5, 4*n-7)] = tt[n-1];
        a[(4*n-5, 4*n-6)] = 1.0;   
    }
    a
}

/// Add to the matrix the coefficients corresponding to the equations of 
/// continuousity of the 1st and 2nd order derivatives.
fn add_cont_equations(a: &mut DMatrix<f64>, i: usize, t: f64, is_linear: bool){
    // linear interpolation
    if is_linear {
        /* So that we justa have a cubic polynomial with a = 0, b = 0 */
        a[(4*i+2, 4*i)]   = 1.0;
        a[(4*i+3, 4*i+1)] = 1.0;
    }
    // cubic spline interpolation
    else {
        let t_squared = t * t;
        /* Slope continuousity */
        a[(4*i+2, 4*i)]   = 3.0*t_squared;
        a[(4*i+2, 4*i+1)] = 2.0*t;
        a[(4*i+2, 4*i+2)] = 1.0;

        a[(4*i+2, 4*i+4)] = -3.0*t_squared;
        a[(4*i+2, 4*i+5)] = -2.0*t;
        a[(4*i+2, 4*i+6)] = -1.0;
    
        /* Concavity continousity */
        a[(4*i+3, 4*i)]   = 6.0*t;
        a[(4*i+3, 4*i+1)] = 2.0;

        a[(4*i+3, 4*i+4)] = -6.0*t;
        a[(4*i+3, 4*i+5)] = -2.0;
    }
}

/// Interpolates points (t[i], x[i]) with a cubic
/// spline and returns it. 
/// May fail if some timestamps are equal, thus trying to
/// inverse a matrix which determinant is 0.
#[allow(dead_code)]
pub fn interpolate (xx: &Vec<f64>, tt: &Vec<f64>, mm: &Vec<bool>) -> Spline {
    assert_eq!(xx.len(), tt.len());
    let n = xx.len();
    let mut s = Spline {parts: Vec::new(),
                        changes: Vec::new(),
                        current: 0, start: xx[0],
                        end: xx[xx.len()-1]};

    let a: DMatrix<f64> = matrix_for(&tt, &mm);
    let mut b: DVector<f64> = DVector::zeros(4*(n-1));

    // Building `b` vector:
    for i in 0..(n-1) {
        s.changes.push(tt[i]);

        b[4*i]   = xx[i];
        b[(4*i+1)] = xx[(i+1)];
    };
    s.changes.push(tt[n-1]);

    // println!("{}{}", a, b);
    let dec = a.lu();// critical point lol
    let x = dec.solve(&b).expect("Computation of spline's coefficients failed !");

    // println!("{}", x);
    for i in 0..(n-1) {
        s.parts.push(SplinePart{
            a: x[4*i],
            b: x[4*i+1],
            c: x[4*i+2],
            d: x[4*i+3],
        });
    }
    s
}

/// Iterpolates sets of points (t[i], x[i][j]) (where i is constant),
/// that share the same timestamps in t. In particular, the coordinates
/// of points to interpolate for the drawing share the same t.
#[allow(dead_code)]
pub fn interpolate_coords(xxx: Vec<Vec<f64>>, tt: &Vec<f64>, mm: &Vec<bool>)
 -> Vec<Spline> {
    let n = tt.len();
    let count = xxx.len();
    let a = matrix_for(&tt, &mm);
    // println!("A = {}", a);
    let dec = a.lu();

    let mut ss = Vec::new();

    for i in 0..count {
        assert_eq!(n, xxx[i].len());
        ss.push(  Spline {parts: Vec::new(),
                          changes: tt.clone(),
                          current: 0,
                          start: tt[0],  end: tt[n-1]    });
    
        let mut b: DVector<f64> = DVector::zeros(4*(n-1));
        // Building `b` vector:
        for j in 0..(n-1) {
            b[4*j]   = xxx[i][j];
            b[(4*j+1)] = xxx[i][(j+1)];
        };
        // println!("xx: {:?} b: {}", xxx[i], b);

        let x = dec.solve(&b).expect("Computation of spline's coefficients failed !");
        // println!("(i={}) sol = {}", i, x);  
        for j in 0..(n-1) {
            ss[i].parts.push(SplinePart{
                a: x[4*j],
                b: x[4*j+1],
                c: x[4*j+2],
                d: x[4*j+3],
            });
        }
    };
    ss
}

/// Evaluates the spline at position x.
#[allow(dead_code)]
pub fn eval(spline: & Spline, x: f64) -> f64 {
    if x < spline.start || x > spline.end {return 0 as f64}
    let mut npart = 0;
    for c in spline.changes() {
        if x < c {break;}
        npart += 1;
    }
    eval_part(spline.parts[npart-1], x)
}

/// Compute value of cubic polynomial hold in the splinepart
/// at position x.
fn eval_part(part: SplinePart, x: f64) -> f64 {
    part.a * x.powi(3) + part.b * x.powi(2) + part.c * x + part.d
}

impl Iterator for Spline {
    type Item = SplinePart;

    fn next (&mut self) -> Option<SplinePart> {
        if self.current < self.parts.len() {
            self.current += 1;
            Some(self.parts[self.current-1])
        }
        else {
            self.current = 0;
            None
        }
    }
}

impl fmt::Display for Spline {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut i: usize = 0;
        for sp in self.parts.clone() {
            write!(f, "Part {}) [t: {:+.4e}]  {:+.4e}x³ {:+.4e}x² {:+.4e}x {:+.4e}\n", 
                        i+1, self.changes[i], sp.a, sp.b, sp.c, sp.d)?;
            i += 1;
        };
        Ok(())
    }
}

#[allow(dead_code)]
impl Spline {
    pub fn start(&self) -> f64 {self.start}
    pub fn end(&self)   -> f64 {self.end}
    pub fn part(&self, i: usize) -> SplinePart {self.parts[i]}
    pub fn changes(&self) -> Vec<f64> {self.changes.clone()}
    pub fn num_parts(&self) -> usize {self.parts.len()}
    pub fn eval(&self, t: f64) -> f64 {eval(self, t)}
}