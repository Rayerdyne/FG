extern crate nalgebra as na;

use na::{DMatrix, DVector};

pub fn hello_spline () {
    println!("I'm in spline !");
}

#[allow(dead_code)]
pub struct Spline {
    parts: Vec<SplinePart>,
    changes: Vec<f64>,
    current: usize,
    start: f64,
    end: f64,
}

#[allow(dead_code)]
#[derive(Copy, Clone)]
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

    let mut a: DMatrix<f64> = DMatrix::zeros(4*n, 4*n);
    let mut b: DVector<f64> = DVector::zeros(4*n);

    let mut t_squared: f64;
    let mut t_cubed: f64;

    // Building `a` matrix and `b` vector:
    for i in 0..n {
        s.changes[i] = tt[i];

        b[4*i]   = xx[i];
        b[4*i+1] = xx[i+1];

        t_squared = tt[i] * tt[i];
        t_cubed = t_squared * tt[i];

        a[(4*i, 4*i)]   = t_cubed;
        a[(4*i, 4*i+1)] = t_squared;
        a[(4*i, 4*i+2)] = tt[i];
        a[(4*i, 4*i+3)] = 1 as f64;

        a[(4*i+1, 4*i)]   = t_cubed;
        a[(4*i+1, 4*i+1)] = t_squared;
        a[(4*i+1, 4*i+2)] = tt[i];
        a[(4*i+1, 4*i+3)] = 1 as f64;

        a[(4*i+2, 4*i)] =   (3 as f64)*t_squared;
        a[(4*i+2, 4*i+1)] = (2 as f64)*tt[i];
        a[(4*i+2, 4*i+2)] = 1 as f64;

        a[(4*i+3, 4*i)] =   (3 as f64)*t_squared;
        a[(4*i+3, 4*i+1)] = (2 as f64)*tt[i];
        a[(4*i+3, 4*i+2)] = 1 as f64;

        if i > 0 {
            a[(4*i+2, 4*i-4)] = -(3 as f64)*t_squared;
            a[(4*i+2, 4*i-3)] = -(2 as f64)*tt[i];
            a[(4*i+2, 4*i-2)] = -1 as f64;
        }

        if i < n-1 {
            a[(4*i+3, 4*i-4)] = -(3 as f64)*t_squared;
            a[(4*i+3, 4*i-3)] = -(2 as f64)*tt[i];
            a[(4*i+3, 4*i-2)] = -1 as f64;
        }
    }

    let dec = a.lu();// critical point lol
    let x = dec.solve(&b).expect("Computation of spline's coefficients failed !");

    for i in 0..n {
        s.parts[i].a = x[4*i];
        s.parts[i].b = x[4*i+1];
        s.parts[i].c = x[4*i+2];
        s.parts[i].d = x[4*i+3];
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