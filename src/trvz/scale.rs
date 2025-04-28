use super::align::{Align, SegType};
use super::params::Color;
use pipeplot::{Pipe, Seg, Shape};

pub fn get_scale(mut xpos: u32, ypos: u32, height: u32, align: &Align) -> Pipe {
    let lf_len = align
        .iter()
        .filter(|op| op.seg_type == SegType::LeftFlank)
        .map(|op| op.width)
        .sum::<usize>();
    xpos += lf_len as u32;

    let allele_len = align
        .iter()
        .filter(|op| op.seg_type != SegType::LeftFlank && op.seg_type != SegType::RightFlank)
        .map(|op| op.width)
        .sum::<usize>();

    let mut label = allele_len.to_string();
    label += "bp";
    let segs = vec![Seg {
        width: allele_len as u32,
        color: Color::Black.to_string(),
        shape: Shape::DoubleArrow(Some(label)),
    }];

    Pipe {
        xpos,
        ypos,
        height,
        segs,
        bands: Vec::new(),
        outline: false,
    }
}
