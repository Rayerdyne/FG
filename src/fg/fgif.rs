extern crate gif;

use gif::{Frame, Encoder, Repeat, SetParameter};
use std::fs::File;
use std::io::Error;
use std::f64::{self, consts::PI};
use std::borrow::Cow;

use super::fourier::CoeffsSet;

#[allow(dead_code)]
struct MyGif<'a> {
    encoder: Encoder<&'a mut File>,
    width: u16,
    height: u16,
}

#[allow(dead_code)]
impl MyGif<'_> {
    fn new<'a> (file: &'a mut File, w: u16, h: u16, global_palette: &[u8]) -> MyGif<'a> {
        let mut encoder = Encoder::new(file, w, h, global_palette).unwrap();
        encoder.set(Repeat::Infinite).unwrap();

        MyGif {
            width: w,           
            height: h,
            encoder: encoder,
        }
    }

    fn write_frame(&mut self, t: &[u8]) {
        let mut frame = Frame::default();
        frame.width  = self.width;
        frame.height = self.height;
        frame.buffer = Cow::Borrowed(&*t);
        self.encoder.write_frame(&frame).unwrap();
    }
}

#[allow(dead_code)]
fn draw_line(xi: usize, yi: usize, xf: usize, yf: usize, color: u8,
             tab: &mut [u8], tabw: usize, tabh: usize) {
    
    let xi2 = if xi >= tabw {  tabw-1  } else {  xi  };
    let xf2 = if xf >= tabw {  tabw-1  } else {  xf  };
    let yi2 = if yi >= tabh {  tabh-1  } else {  yi  };
    let yf2 = if yf >= tabh {  tabh-1  } else {  yf  };

    assert!(xi2 < tabw && xf2 < tabw);
    assert!(yi2 < tabh && yf2 < tabh);
    let (x1, x2, x_inversed) = if xf2 < xi2 { (xf2, xi2, true)   } 
                               else         { (xi2, xf2, false)  };
    let (y1, y2, y_inversed) = if yf2 < yi2 { (yf2, yi2, true)   }
                               else         { (yi2, yf2, false)  };
    
    if x1 == x2 {
        if y1 == y2 {   tab[y1*tabw + x1] = color; 
                        return }
        for i in y1..=y2 {   tab[i*tabw + x1] = color;   }
        return
    }
    if y1 == y2 {
        for i in x1..=x2 {   tab[y1*tabw + i] = color;   }
        return

    }

    let mut y = y1;
    let (dx, dy) = ((x2 as i32 - x1 as i32), (y2 as i32 - y1 as i32));
    let mut e: i32 = -dx;
    let (ex, ey) = (2*dy, -2*dx);


    if !x_inversed {
        if !y_inversed { // quadrant 4
            for x in x1..x2 {
                tab[y*tabw + x] = color;

                e += ex;
                while e >= 0 {
                    y += 1;
                    e += ey;
                    if e >= 0 {   tab[y*tabw + x] = color;   }
                }
            }
        }
        else { // quadrant 1 
            for x in x1..x2 {
                tab[(y1+y2-y)*tabw + x] = color;

                e += ex;
                while e >= 0 {
                    y += 1;
                    e += ey;
                    if e >= 0 {   tab[(y1+y2-y)*tabw + x] = color;   }
                }
            }
        }
    }
    else { 
        if !y_inversed { // quadrant 3
            for x in x1..x2 {
                tab[y*tabw + x1+x2-x] = color;

                e += ex;
                while e >= 0 {
                    y += 1;
                    e += ey;
                    if e >= 0 {   tab[y*tabw + x1+x2-x] = color;   }
                }
            }
        }
        else { // quadrant 2
            for x in x1..x2 {
                tab[(y1+y2-y)*tabw + x1+x2-x] = color;

                e += ex;
                while e >= 0 {
                    y += 1;
                    e += ey;
                    if e >= 0 {   tab[(y1+y2-y)*tabw + x1+x2-x] = color;   }
                }
            }
        }
    }

}

/* see https://fr.wikipedia.org/wiki/Algorithme_de_trac%C3%A9_de_segment_de_Bresenham */

pub fn draw_fourier_coeff(coeffs: CoeffsSet, filename: &str,
                 w: usize, h: usize, global_palette: &[u8]) -> Result<(), Error> {
    gotest(20, 20, 200, 200);
    let n = coeffs.ppos.len();
    assert_eq!(n, coeffs.nneg.len());

    let mut output = File::create(filename)?;
    let mut gif = MyGif::new(&mut output, w as u16, h as u16, global_palette);
    
    let vect = vec![0; (w*h) as usize];
    let mut tab_drawing: Box<[u8]> = vect.into_boxed_slice();
    println!("w: {}, h: {}, {}", w as u16, h as u16, (w*h) as usize);
    
    let mut t: f64 = 0.0_f64;
    let max = 2.0_f64*PI;
    while t < max {
        // keep what's already drawed
        let mut tab_lines = tab_drawing.clone();
        
        let mut x1: f64 = (w as f64) / 2.0_f64;
        let mut y1: f64 = (h as f64) / 2.0_f64;
        let mut x2: f64 = 0.0_f64;
        let mut y2: f64 = 0.0_f64;

        let mut k_f64: f64 = 1.0_f64;
        for k in 0..n {
            let coeff_pn = (coeffs.ppos[k], coeffs.nneg[k]);
            for (c, neg) in vec![(coeff_pn.0, false),
                                 (coeff_pn.1, true) ] {
                let sin1 = if neg {  -(k_f64*t).sin()  }
                           else   {   (k_f64*t).sin()  };
                let cos1 = (k_f64*t).cos();
                //   (a+ib)*(cos + i sin)
                // = a cos - b sin + i (a sin + b cos)
                x2 = x1 + (c.a*cos1 - c.b*sin1);
                y2 = y1 - (c.a*sin1 + c.b*cos1);
                // Y axis is multiplied by -1 to make the circle drawed anticlockwise 
                draw_line(
                    if x1 < 0.0_f64 { 0 } else { x1 as usize },
                    if y1 < 0.0_f64 { 0 } else { y1 as usize },
                    if x2 < 0.0_f64 { 0 } else { x2 as usize }, 
                    if y2 < 0.0_f64 { 0 } else { y2 as usize },
                    1, &mut *tab_lines, w, h);
                x1 = x2;
                y1 = y2;
            }

            k_f64 += 1.0_f64;
        }
        let (xx, yy) = (if x2 < 0.0_f64 { 0 } else { x2 as usize },
                        if y2 < 0.0_f64 { 0 } else { y2 as usize } );
        draw_line(xx, yy, xx, yy, 1, 
                    &mut *tab_drawing, w, h);

        gif.write_frame(&mut *tab_lines);
        t += 0.05;
    };

    Ok(())
}

#[allow(dead_code)]
pub fn gotest(a: usize, b: usize, c: usize, d: usize) {
    let mut f = File::create("hello.gif").expect("couldn't create file");

    let (w, h) = (300 as usize, 200 as usize);
    let global_palette = &[0x00, 0x00, 0xFF,
                           0xFF, 0x00, 0x00];
    let mut gif = MyGif::new(&mut f, w as u16, h as u16, global_palette);

    let vect = vec![0; w*h];
    let mut tab: Box<[u8]> = vect.into_boxed_slice();
    draw_line(a, b, c, d, 1, &mut *tab, w, h);
    draw_line(100, 100, 150, 0, 1, &mut *tab, w, h);
    draw_line(100, 100, 50, 200, 1, &mut *tab, w, h);
    gif.write_frame(&*tab);
}