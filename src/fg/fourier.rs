use super::spline::*;
use std::f64::{self, consts::PI};

#[allow(dead_code)]
pub fn compute_fourier_coeff(sx: Spline, sy: Spline, n: usize) -> (Vec<[f64; 2]>, Vec<[f64; 2]>) {
    let (t_i, t_f) = (sx.start(), sx.end());
    assert_eq!(t_i, sy.start());
    assert_eq!(t_f, sy.end());
    assert_eq!(sx.num_parts(), sy.num_parts());
    let period = t_f - t_i;

    let omega_0 = 2.0_f64*PI / period;
    
    let mut pp: Vec<[f64; 2]> = vec![[0 as f64, 0 as f64]; n];
    let mut nn: Vec<[f64; 2]> = vec![[0 as f64, 0 as f64]; n];

    let changes = sx.changes();
    let part_count = changes.len()-1;

    let mut vx = CubicIntegrator::new(t_i, omega_0);
    let mut vy = CubicIntegrator::new(t_i, omega_0);

    for i in 0..part_count {
        let t2 = changes[i+1];
        vx.next_step(sx.part(i), t2);
        vy.next_step(sy.part(i), t2);
        add_splines_contributions(&mut pp, &mut nn, &vx, &vy);
    };

    (pp, nn)
}

#[allow(dead_code)]
fn add_splines_contributions(pp: &mut Vec<[f64; 2]>, nn: &mut Vec<[f64; 2]>, 
                                  vx: &CubicIntegrator, vy: &CubicIntegrator) {
    let n = pp.len();
    assert_eq!(n, nn.len());

    for i in 1..n {
        //contribution of X spline
        let p_contr_x = integral_12(vx, i, false);
        let n_contr_x = integral_12(vx, i, true);

        pp[i][0] += p_contr_x[0];         nn[i][0] += n_contr_x[0];
        pp[i][1] += p_contr_x[1];         nn[i][1] += n_contr_x[1];

        //contribution of Y spline
        let p_contr_y = integral_12(vy, i, false);
        let n_contr_y = integral_12(vy, i, true);

        pp[i][0] += -p_contr_y[1];        nn[i][0] += -n_contr_y[1];
        pp[i][1] += p_contr_y[0];         nn[i][1] += n_contr_y[0];
        //    \>   i*(a+ib) = -b+ia
        }
    // ok this is not optimal.
    // but it will be ok
}

fn integral_12(v: &CubicIntegrator, k_index: usize, negative: bool) -> [f64; 2] {
    let k: f64 = if negative { (k_index as f64)*(-1.0_f64)}
                 else {k_index as f64};
    let k_sq = k.powf(2.0_f64);
    let k_cu = k.powf(3.0_f64);
    let k_fo = k.powf(4.0_f64);
    let arg1 = -v.omega_0*k*v.t1;
    let arg2 = -v.omega_0*k*v.t2;

    let cos1 = arg1.cos();       let sin1 = arg1.sin();
    let cos2 = arg2.cos();       let sin2 = arg2.sin();

    let term_1_re = (-sin1 * v.r1_1 / k ) + 
                    ( cos1 * v.r1_2 / k_sq) +
                    ( sin1 * v.r1_2 / k_cu) +
                    (-cos1 * v.r1_2 / k_fo);
                    
    let term_1_im = ( cos1 * v.r1_1 / k ) + 
                    ( sin1 * v.r1_2 / k_sq) +
                    (-cos1 * v.r1_2 / k_cu) +
                    (-sin1 * v.r1_2 / k_fo);

    let term_2_re = (-sin2 * v.r2_1 / k ) + 
                    ( cos2 * v.r2_2 / k_sq) +
                    ( sin2 * v.r2_2 / k_cu) +
                    (-cos2 * v.r2_2 / k_fo);

    let term_2_im = ( cos2 * v.r2_1 / k ) + 
                    ( sin2 * v.r2_2 / k_sq) +
                    (-cos2 * v.r2_2 / k_cu) +
                    (-sin2 * v.r2_2 / k_fo);
    
    [term_2_re - term_1_re, term_2_im - term_1_im]
}

struct CubicIntegrator {
    r1_1: f64,      r2_1: f64,
    r1_2: f64,      r2_2: f64,
    r1_3: f64,      r2_3: f64,
    r1_4: f64,      r2_4: f64,

    t1: f64,        t2: f64,
    t1_sq: f64,     t2_sq: f64,
    t1_cu: f64,     t2_cu: f64,
    a: f64,         b: f64,         c: f64,      d: f64,

    omega_0: f64,   omega_0_sq: f64,    omega_0_cu: f64,    omega_0_fo: f64
}

impl CubicIntegrator {
    fn next_step(&mut self, sp: SplinePart, t2: f64) {  
        self.a = sp.geta();           self.b = sp.getb();
        self.c = sp.getc();           self.d = sp.getd();

        self.t1 = self.t2;
        self.t1_sq = self.t2_sq;
        self.t1_cu = self.t2_cu;
        self.t2 = t2;
        self.t2_sq = t2.powf(2.0_f64);
        self.t2_cu = t2.powf(3.0_f64);

        let c1_1 = self.a*self.t1_cu +     self.b*self.t1_sq +     self.c*self.t1 + self.d;
        let c1_2 =                     3.0*self.a*self.t1_sq + 2.0*self.b*self.t1 + self.c;
        let c1_3 =                                             6.0*self.a*self.t1 + 2.0*self.b;
        let c1_4 =                                                                  6.0*self.a;
        
        let c2_1 = self.a*self.t2_cu +     self.b*self.t2_sq +     self.c*self.t2 + self.d;
        let c2_2 =                     3.0*self.a*self.t2_sq + 2.0*self.b*self.t2 + self.c;
        let c2_3 =                                             6.0*self.a*self.t2 + 2.0*self.b;
        let c2_4 =                                                                  6.0*self.a;


        self.r1_1 = c1_1 / self.omega_0;                self.r2_1 = c2_1 / self.omega_0;    
        self.r1_2 = c1_2 / self.omega_0_sq;             self.r2_2 = c2_2 / self.omega_0_sq;
        self.r1_3 = c1_3 / self.omega_0_cu;             self.r2_3 = c2_3 / self.omega_0_cu;
        self.r1_4 = c1_4 / self.omega_0_fo;             self.r2_4 = c2_4 / self.omega_0_fo;  
    }

    fn new(t_i: f64, omega_0: f64) -> CubicIntegrator {
        CubicIntegrator {
            a: 0.0_f64,              b: 0.0_f64,
            c: 0.0_f64,              d: 0.0_f64,

            t1: 0.0_f64,             t2: t_i,
            t1_sq: 0.0_f64,          t2_sq: t_i.powf(2.0_f64),
            t1_cu: 0.0_f64,          t2_cu: t_i.powf(3.0_f64),

            r1_1: 0.0_f64,          r2_1: 0.0_f64,
            r1_2: 0.0_f64,          r2_2: 0.0_f64,
            r1_3: 0.0_f64,          r2_3: 0.0_f64,
            r1_4: 0.0_f64,          r2_4: 0.0_f64,
            omega_0: omega_0,                   omega_0_sq: omega_0.powf(2.0_f64),
            omega_0_cu: omega_0.powf(3.0_f64),  omega_0_fo: omega_0.powf(4.0_f64)  
        }
    }
}
