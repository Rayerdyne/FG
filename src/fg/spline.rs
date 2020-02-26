extern crate nalgebra as na;

use na::{DMatrix, DVector};
use std::fmt;

#[allow(dead_code)]
#[derive(Clone, Debug)]
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

fn matrix_for (tt: &Vec<f64>) -> DMatrix<f64> {
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
        a[(4*i, 4*i+3)] = 1.0_f64;

        a[(4*i+1, 4*i)]   = t2_cubed;
        a[(4*i+1, 4*i+1)] = t2_squared;
        a[(4*i+1, 4*i+2)] = tt[i+1];
        a[(4*i+1, 4*i+3)] = 1.0_f64;
    
        if i < n-2 {
            /* Slope continuousity */
            a[(4*i+2, 4*i)]   = 3.0_f64*t2_squared;
            a[(4*i+2, 4*i+1)] = 2.0_f64*tt[i+1];
            a[(4*i+2, 4*i+2)] = 1.0_f64;

            a[(4*i+2, 4*i+4)] = -3.0_f64*t2_squared;
            a[(4*i+2, 4*i+5)] = -2.0_f64*tt[i+1];
            a[(4*i+2, 4*i+6)] = -1.0_f64;
        
            /* Concavity continousity */
            a[(4*i+3, 4*i)]   = 6.0_f64*tt[i+1];
            a[(4*i+3, 4*i+1)] = 2.0_f64;

            a[(4*i+3, 4*i+4)] = -6.0_f64*tt[i+1];
            a[(4*i+3, 4*i+5)] = -2.0_f64;
        }

        /*if i > 0 {
            a[(4*i+2, 4*i-4)] = -3.0_f64*t1_squared;
            a[(4*i+2, 4*i-3)] = -2.0_f64*tt[i];
            a[(4*i+2, 4*i-2)] = -1.0_f64;
        }

        if i < n-2 {
            a[(4*i+3, 4*i+4)] = -3.0_f64*t2_squared;
            a[(4*i+3, 4*i+5)] = -2.0_f64*tt[i+1];
            a[(4*i+3, 4*i+6)] = -1.0_f64;
        }*/
    }

    /* Starting & ending slopes are 0 */
    let ti_squared = tt[0] * tt[0];
    a[(4*n-6, 0)] = ti_squared;
    a[(4*n-6, 1)] = tt[0];
    a[(4*n-6, 2)] = 1.0_f64;
    
    let tf_squared = tt[n-1] * tt[n-1];
    a[(4*n-5, 4*n-8)] = tf_squared;
    a[(4*n-5, 4*n-7)] = tt[n-1];
    a[(4*n-5, 4*n-6)] = 1.0_f64;   
    a
}

#[allow(dead_code)]
pub fn interpolate (xx: Vec<f64>, tt: Vec<f64>) -> Spline {
    assert_eq!(xx.len(), tt.len());
    let n = xx.len();
    let mut s = Spline {parts: Vec::new(),
                        changes: Vec::new(),
                        current: 0, start: xx[0],
                        end: xx[xx.len()-1]};

    let a: DMatrix<f64> = matrix_for(&tt);
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

#[allow(dead_code)]
pub fn interpolate_coords(xxx: Vec<Vec<f64>>, tt: Vec<f64>) -> Vec<Spline> {
    let n = tt.len();
    let count = xxx.len();
    let a = matrix_for(&tt);
    println!("A = {}", a);
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
        println!("xx: {:?} b: {}", xxx[i], b);

        let x = dec.solve(&b).expect("Computation of spline's coefficients failed !");
        println!("(i={}) sol = {}", i, x);  
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
}

#[allow(dead_code)]
impl SplinePart {
    pub fn geta(&self) -> f64 {self.a}
    pub fn getb(&self) -> f64 {self.b}
    pub fn getc(&self) -> f64 {self.c}
    pub fn getd(&self) -> f64 {self.d}
}