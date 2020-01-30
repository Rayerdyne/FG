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

    fn write_frame(t: &[u8], n: usize, gif: &mut MyGif) {
        assert_eq!(n, (gif.width * gif.height) as usize);
        let mut frame = Frame::default();
        frame.width  = gif.width;
        frame.height = gif.height;
        frame.buffer = Cow::Borrowed(&*t);
        gif.encoder.write_frame(&frame).unwrap();
    }
}
