use super::align::SegType;
use std::{collections::HashMap, fmt};

pub type ColorMap = HashMap<SegType, Color>;

pub fn pick_colors(motifs: &[String]) -> ColorMap {
    let tr_colors = [
        Color::Blue,
        Color::Purple,
        Color::Orange,
        Color::Pink,
        Color::Yellow,
        Color::Green,
        Color::Red,
        Color::Khaki,
        Color::PaleRed,
        Color::PaleBlue,
    ];

    let mut colors = HashMap::new();
    colors.insert(SegType::LeftFlank, Color::Teal);
    colors.insert(SegType::RightFlank, Color::Teal);

    for index in 0..motifs.len() {
        let color = tr_colors[index % tr_colors.len()].clone();
        colors.insert(SegType::Tr(index), color);
    }
    colors.insert(SegType::Tr(motifs.len()), Color::Gray);

    colors
}

pub fn get_meth_colors(motifs: &[String]) -> ColorMap {
    let mut colors = HashMap::new();
    colors.insert(SegType::LeftFlank, Color::Teal);
    colors.insert(SegType::RightFlank, Color::Teal);

    for index in 0..motifs.len() + 1 {
        colors.insert(SegType::Tr(index), Color::Gray);
    }

    colors
}

#[derive(Debug, PartialEq, Clone)]
pub enum Color {
    Purple,
    Blue,
    Orange,
    Teal,
    Gray,
    LightGray,
    Black,
    Green,
    Pink,
    Yellow,
    Red,
    Khaki,
    PaleRed,
    PaleBlue,
    Grad(f64),
}

impl fmt::Display for Color {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Color::Purple => write!(formatter, "#814ED1"),
            Color::Blue => write!(formatter, "#1383C6"),
            Color::Orange => write!(formatter, "#E16A2C"),
            Color::Teal => write!(formatter, "#009CA2"),
            Color::Gray => write!(formatter, "#BABABA"),
            Color::LightGray => write!(formatter, "#D1D1D1"),
            Color::Black => write!(formatter, "#000000"),
            Color::Pink => write!(formatter, "#ED3981"),
            Color::Yellow => write!(formatter, "#EFCD17"),
            Color::Green => write!(formatter, "#009D4E"),
            Color::Red => write!(formatter, "#E3371E"),
            Color::Khaki => write!(formatter, "#F0E68C"),
            Color::PaleRed => write!(formatter, "#FF4858"),
            Color::PaleBlue => write!(formatter, "#46B2E8"),
            Color::Grad(value) => write!(formatter, "{}", get_gradient(*value)),
        }
    }
}

fn get_gradient(value: f64) -> String {
    let blue: (u8, u8, u8) = (0, 73, 255);
    let red: (u8, u8, u8) = (255, 0, 0);
    let mix_red = (blue.0 as f64 * (1.0 - value) + red.0 as f64 * value).round() as u8;
    let mix_green = (blue.1 as f64 * (1.0 - value) + red.1 as f64 * value).round() as u8;
    let mix_blue = (blue.2 as f64 * (1.0 - value) + red.2 as f64 * value).round() as u8;

    format!("#{:02X}{:02X}{:02X}", mix_red, mix_green, mix_blue)
}
