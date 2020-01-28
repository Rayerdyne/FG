extern crate nalgebra as na;

use na::{DMatrix, DVector};
use std::fmt;

#[allow(dead_code)]
#[derive(Debug)]
pub struct Spline {
    parts: Vec<SplinePart>,
    changes: Vec<f64>,
    current: usize,
    start: f64,
    end: f64,
}

#[allow(dead_code)]
#[derive(Copy, Clone, Debug)]
pub struct SplinePart {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
}

#[allow(dead_code)]
pub fn interpolate (xx: Vec<f64>, tt: Vec<f64>) -> Spline {
    assert_eq!(xx.len(), tt.len());
    let n = xx.len();
    let mut s = Spline {parts: Vec::new(),
                        changes: Vec::new(),
                        current: 0, start: xx[0],
                        end: xx[xx.len()-1]};

    let mut a: DMatrix<f64> = DMatrix::zeros(4*(n-1), 4*(n-1));
    let mut b: DVector<f64> = DVector::zeros(4*(n-1));

    let mut t1_squared: f64 = tt[0] * tt[0];
    let mut t1_cubed:   f64 = t1_squared * tt[0];
    let mut t2_squared: f64 = tt[1] * tt[1];
    let mut t2_cubed:   f64 = t2_squared * tt[1];

    // Building `a` matrix and `b` vector:
    for i in 0..(n-1) {
        s.changes.push(tt[i]);

        b[4*i]   = xx[i];
        b[(4*i+1)] = xx[(i+1)];

        if i >= 1 {
            t1_squared = t2_squared;
            t1_cubed = t2_cubed;
        }
        t2_squared = tt[i+1] * tt[i+1];
        t2_cubed = t2_squared * tt[i+1];

        a[(4*i, 4*i)]   = t1_cubed;
        a[(4*i, 4*i+1)] = t1_squared;
        a[(4*i, 4*i+2)] = tt[i];
        a[(4*i, 4*i+3)] = 1 as f64;

        a[(4*i+1, 4*i)]   = t2_cubed;
        a[(4*i+1, 4*i+1)] = t2_squared;
        a[(4*i+1, 4*i+2)] = tt[i+1];
        a[(4*i+1, 4*i+3)] = 1 as f64;

        a[(4*i+2, 4*i)] =   (3 as f64)*t1_squared;
        a[(4*i+2, 4*i+1)] = (2 as f64)*tt[i];
        a[(4*i+2, 4*i+2)] = 1 as f64;

        a[(4*i+3, 4*i)] =   (3 as f64)*t2_squared;
        a[(4*i+3, 4*i+1)] = (2 as f64)*tt[i+1];
        a[(4*i+3, 4*i+2)] = 1 as f64;

        if i > 0 {
            a[(4*i+2, 4*i-4)] = -(3 as f64)*t1_squared;
            a[(4*i+2, 4*i-3)] = -(2 as f64)*tt[i];
            a[(4*i+2, 4*i-2)] = -1 as f64;
        }

        if i < n-2 {
            a[(4*i+3, 4*i+4)] = -(3 as f64)*t2_squared;
            a[(4*i+3, 4*i+6)] = -(2 as f64)*tt[i+1];
            a[(4*i+3, 4*i+7)] = -1 as f64;
        }
    }

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

#[allow(dead_code)]
pub fn eval(spline: Spline, x: f64) -> f64 {
    let mut npart = 0;
    for c in spline.changes {
        npart += 1;
        if c < x {break;}
    }
    if x < spline.start || x < spline.end          {0 as f64}
    else {eval_part(spline.parts[npart], x)}
}

#[allow(dead_code)]
pub fn eval_part(part: SplinePart, x: f64) -> f64 {
    part.a * x * x * x + part.b * x * x + part.c * x + part.d
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
            write!(f, "Part {}) [t: {}]  {}x³ + {}x² + {}x + {}", 
                        i+1, self.changes[i], sp.a, sp.b, sp.c, sp.d)?;
            i += 1;
        };
        Ok(())
    }
}