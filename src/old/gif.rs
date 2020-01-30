pub struct FGif {
    header: FGifHeader,
    body: Vec<FGifImage>,
}

struct FGifHeader {
    width: u16,         // of dilplay screen
    height: u16,
    packed: u8,         // screen & color map information
    background: u8,     // its index
    aspect_ratio: u8,   // pixel aspect ratio
}

struct FGifImage {
    pos_x: u16,
    pos_y: u16,
    width: u16,
    height: u16,
    packed: u8,
}

pub struct ColorTable {
    r: u8,
    g: u8,
    b: u8,
}

pub enum WriteError {

}

impl FGif {
    fn write_to (filename: &str) -> Result<(), WriteError> {
        let mut v = Vec::<u8>::new();
        //constant header
        v.extend(b"GIF89a".iter().cloned());

        Ok(())
    }
}

/* 
        FGifHeader::packed 
Bits 0-2 	Size of the Global Color Table
Bit 3 	Color Table Sort Flag
Bits 4-6 	Color Resolution
Bit 7 	Global Color Table Flag (y-n)

see: https://www.fileformat.info/format/gif/egff.htm
*/

/*
        FGifImage::packed
Bit 0 	Local Color Table Flag (y-n)
Bit 1 	Interlace Flag
Bit 2 	Sort Flag
Bits 3-4 	Reserved
Bits 5-7 	Size of Local Color Table Entry
*/