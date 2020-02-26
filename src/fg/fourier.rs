use super::spline::*;
use std::f64::{self, consts::PI};

#[derive(Debug, Clone, Copy)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

#[derive(Debug)]
pub struct CoeffsSet {
    pub ppos: Vec<Complex>,
    pub nneg: Vec<Complex>,
    //doubled because vector
}

impl std::ops::AddAssign for Complex {
    fn add_assign(&mut self, c: Complex) {
        self.re += c.re;
        self.im += c.im;
    }
}

impl Complex {
    fn times_i(&mut self) -> Complex {
        let temp = self.re;
        self.re = -self.im;
        self.im = temp;
        *self
    }
}

impl CoeffsSet {
    fn new(n: usize) -> CoeffsSet {
        CoeffsSet {
            ppos: Vec::with_capacity(n),
            nneg: Vec::with_capacity(n),
        }
    }
}

#[allow(dead_code)]
pub fn compute_fourier_coeff(sx: Spline, sy: Spline, n: usize) -> CoeffsSet {
    let (t_i, t_f) = (sx.start(), sx.end());
    assert_eq!(t_i, sy.start());
    assert_eq!(t_f, sy.end());
    assert_eq!(sx.num_parts(), sy.num_parts());

    let period = t_f - t_i;
    let omega_0 = 2.0_f64*PI / period;
    
    let mut coeffs = CoeffsSet::new(n);

    let changes = sx.changes();
    let num_parts = changes.len()-1;

    let mut vx = CubicIntegrator::new(t_i, omega_0);
    let mut vy = CubicIntegrator::new(t_i, omega_0);

    for i in 0..num_parts {
        let t2 = changes[i+1];
        vx.next_step(sx.part(i), t2);
        vy.next_step(sy.part(i), t2);
        add_splines_contributions(&mut coeffs, &vx, &vy);
    };

    coeffs
}

#[allow(dead_code)]
fn add_splines_contributions(coeffs: &mut CoeffsSet, vx: &CubicIntegrator, vy: &CubicIntegrator) {
    let n = coeffs.ppos.len();
    assert_eq!(n, coeffs.nneg.len());

    for i in 1..n {
        //contribution of X spline
        let p_contr_x = integral_12(vx, i, false);
        let n_contr_x = integral_12(vx, i, true);

        coeffs.ppos[i] += p_contr_x;         
        coeffs.nneg[i] += n_contr_x;

        //contribution of Y spline
        let mut p_contr_y = integral_12(vy, i, false);
        let mut n_contr_y = integral_12(vy, i, true);

        coeffs.ppos[i] += p_contr_y.times_i();
        coeffs.nneg[i] += n_contr_y.times_i();
        }
    // ok this is not optimal.
    // but it will be ok
}

fn integral_12(v: &CubicIntegrator, k_index: usize, negative: bool) -> Complex {
    let k: f64 = if negative { (k_index as f64)*(-1.0_f64)}
                 else {k_index as f64};
    let k_sq = k.powf(2.0_f64);
    let k_cu = k.powf(3.0_f64);
    let k_fo = k.powf(4.0_f64);
    let arg1 = -v.omega_0*k*v.t1;
    let arg2 = -v.omega_0*k*v.t2;

    let cos1 = arg1.cos();       let sin1 = arg1.sin();
    let cos2 = arg2.cos();       let sin2 = arg2.sin();

    let term_1_re = (-sin1 * v.r1.m1 / k ) + 
                    ( cos1 * v.r1.m2 / k_sq) +
                    ( sin1 * v.r1.m3 / k_cu) +
                    (-cos1 * v.r1.m4 / k_fo);
                    
    let term_1_im = ( cos1 * v.r1.m1 / k ) + 
                    ( sin1 * v.r1.m2 / k_sq) +
                    (-cos1 * v.r1.m3 / k_cu) +
                    (-sin1 * v.r1.m4 / k_fo);

    let term_2_re = (-sin2 * v.r2.m1 / k ) + 
                    ( cos2 * v.r2.m2 / k_sq) +
                    ( sin2 * v.r2.m3 / k_cu) +
                    (-cos2 * v.r2.m4 / k_fo);

    let term_2_im = ( cos2 * v.r2.m1 / k ) + 
                    ( sin2 * v.r2.m2 / k_sq) +
                    (-cos2 * v.r2.m3 / k_cu) +
                    (-sin2 * v.r2.m4 / k_fo);
    
    Complex {
        re: term_2_re - term_1_re,
        im: term_2_im - term_1_im,
    }
}

#[derive(Clone, Copy)]
struct VarSet {
    m1: f64,
    m2: f64,
    m3: f64,
    m4: f64,
}

impl VarSet {
    fn new_c(t: f64, t_sq: f64, t_cu: f64, sp: SplinePart) -> VarSet {
        VarSet {
            m1: sp.a*t_cu +    sp.b*t_sq +     sp.c*t +     sp.d,
            m2:            3.0*sp.a*t_sq + 2.0*sp.b*t +     sp.c,
            m3:                            6.0*sp.a*t + 2.0*sp.b,
            m4:                                         6.0*sp.a
        }
    }

    fn new_r(omega: f64, omega_sq: f64, omega_cu: f64, omega_fo: f64, c: VarSet) -> VarSet {
        VarSet {
            m1: c.m1 / omega,
            m2: c.m2 / omega_sq,
            m3: c.m3 / omega_cu,
            m4: c.m4 / omega_fo,
        }
    }

    fn new0() -> VarSet {
        VarSet { m1: 0.0_f64, m2: 0.0_f64, m3: 0.0_f64, m4: 0.0_f64}
    }
}

struct CubicIntegrator {
    r1: VarSet,
    r2: VarSet,

    t1: f64,        t2: f64,
    t1_sq: f64,     t2_sq: f64,
    t1_cu: f64,     t2_cu: f64,

    omega_0: f64,   omega_0_sq: f64,    omega_0_cu: f64,    omega_0_fo: f64
}

impl CubicIntegrator {
    fn next_step(&mut self, sp: SplinePart, t2: f64) {
        self.t1 = self.t2;
        self.t1_sq = self.t2_sq;
        self.t1_cu = self.t2_cu;
        self.t2 = t2;
        self.t2_sq = t2.powf(2.0_f64);
        self.t2_cu = t2.powf(3.0_f64);

        let c1 = VarSet::new_c(self.t1, self.t1_sq, self.t1_cu, sp);
        let c2 = VarSet::new_c(self.t2, self.t2_sq, self.t2_cu, sp);

// let c1_1 = self.re*self.t1_cu +     self.im*self.t1_sq +     self.c*self.t1 + self.d;
        // let c1_2 =                     3.0*self.re*self.t1_sq + 2.0*self.im*self.t1 + self.c;
        // let c1_3 =                                             6.0*self.re*self.t1 + 2.0*self.im;
        // let c1_4 =                                                                  6.0*self.re;
        
        // let c2_1 = self.re*self.t2_cu +     self.im*self.t2_sq +     self.c*self.t2 + self.d;
        // let c2_2 =                     3.0*self.re*self.t2_sq + 2.0*self.im*self.t2 + self.c;
        // let c2_3 =                                             6.0*self.re*self.t2 + 2.0*self.im;
        // let c2_4 =                                                                  6.0*self.re;

        self.r1 = VarSet::new_r(self.omega_0, self.omega_0_sq, self.omega_0_cu, self.omega_0_fo, c1);
        self.r2 = VarSet::new_r(self.omega_0, self.omega_0_sq, self.omega_0_cu, self.omega_0_fo, c2);
        // self.r1_1 = c1_1 / self.omega_0;                self.r2_1 = c2_1 / self.omega_0;    
        // self.r1_2 = c1_2 / self.omega_0_sq;             self.r2_2 = c2_2 / self.omega_0_sq;
        // self.r1_3 = c1_3 / self.omega_0_cu;             self.r2_3 = c2_3 / self.omega_0_cu;
        // self.r1_4 = c1_4 / self.omega_0_fo;             self.r2_4 = c2_4 / self.omega_0_fo;  
    }

    fn new(t_i: f64, omega_0: f64) -> CubicIntegrator {
        CubicIntegrator {
            t1: 0.0_f64,             t2: t_i,
            t1_sq: 0.0_f64,          t2_sq: t_i.powf(2.0_f64),
            t1_cu: 0.0_f64,          t2_cu: t_i.powf(3.0_f64),

            r1: VarSet::new0(),
            r2: VarSet::new0(),
            omega_0: omega_0,                   omega_0_sq: omega_0.powf(2.0_f64),
            omega_0_cu: omega_0.powf(3.0_f64),  omega_0_fo: omega_0.powf(4.0_f64)  
        }
    }
}
