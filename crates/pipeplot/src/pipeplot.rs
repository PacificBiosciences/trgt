#[derive(Debug, PartialEq)]
pub enum Shape {
    Rect,
    HLine,
    VLine,
    None,
    Tick(Option<u32>),
    DoubleArrow(Option<String>),
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
pub struct FontConfig {
    pub family: String,
    pub weight: String,
    pub size: String,
}

impl Default for FontConfig {
    fn default() -> Self {
        Self {
            family: "Roboto Mono".to_string(),
            weight: "bold".to_string(),
            size: "14px".to_string(),
        }
    }
}

#[derive(Debug)]
pub struct PipePlot {
    pub pipes: Vec<Pipe>,
    pub legend: Legend,
    pub font: FontConfig,
}

impl PipePlot {
    pub fn set_font_family(&mut self, font_family: &str) {
        self.font.family = font_family.to_owned();
    }
}
