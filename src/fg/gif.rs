extern crate gif;

use gif::{Frame, Encoder, Repeat, SetParameter};
use std::fs::File;
use std::borrow::Cow;

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

    fn write_frame(&mut self, t: &[u8], n: usize) {
        assert_eq!(n, (self.width * self.height) as usize);
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
                 
    assert!(xi < tabh && yi < tabh);
    assert!(yi < tabw && yf < tabw);
    let (x1, x2, x_inversed) = if xf < xi { (xf, xi, true)   } 
                               else {       (xi, xf, false)  };
    let (y1, y2, y_inversed) = if yf < yi { (yf, yi, true)   }
                               else       { (yi, yf, false)  };
    
    if x1 == x2 {
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
            for x in x1..=x2 {
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
            for x in x1..=x2 {
                tab[(yi+yf-y)*tabw + x] = color;

                e += ex;
                while e >= 0 {
                    y += 1;
                    e += ey;
                    if e >= 0 {   tab[y*tabw + x] = color;   }
                }
            }
        }
    }
    else { 
        if !y_inversed { // quadrant 3
            for x in x1..=x2 {
                tab[y*tabw + xi+xf-x] = color;

                e += ex;
                while e >= 0 {
                    y += 1;
                    e += ey;
                    if e >= 0 {   tab[y*tabw + x] = color;   }
                }
            }
        }
        else { // quadrant 2
            for x in x1..=x2 {
                tab[(yi+yf-y)*tabw + xi+xf-x] = color;

                e += ex;
                while e >= 0 {
                    y += 1;
                    e += ey;
                    if e >= 0 {   tab[y*tabw + x] = color;   }
                }
            }
        }
    }

}

/* see https://fr.wikipedia.org/wiki/Algorithme_de_trac%C3%A9_de_segment_de_Bresenham */


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
    gif.write_frame(&*tab, w*h);
}