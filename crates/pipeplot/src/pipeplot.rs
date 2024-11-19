#[derive(Debug, PartialEq)]
pub enum Shape {
    Rect,
    HLine,
    VLine,
    None,
    Tick(Option<u32>),
}

pub type Color = String;

#[derive(Debug, PartialEq)]
pub struct Seg {
    pub width: u32,
    pub color: Color,
    pub shape: Shape,
}

#[derive(Debug)]
pub struct Band {
    pub pos: u32, // Position relative to pipe's start
    pub width: u32,
    pub color: Color,
}

#[derive(Debug)]
pub struct Pipe {
    pub xpos: u32,
    pub ypos: u32,
    pub height: u32,
    pub segs: Vec<Seg>,
    pub bands: Vec<Band>,
    pub outline: bool,
}

#[derive(Debug)]
pub struct Legend {
    pub xpos: u32,
    pub ypos: u32,
    pub height: u32,
    pub labels: Vec<(String, String)>,
}

#[derive(Debug)]
pub struct PipePlot {
    pub pipes: Vec<Pipe>,
    pub legend: Legend,
}
