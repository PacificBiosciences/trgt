use crate::trvz::{
    align::AlignInfo,
    locus::{BaseLabel, Locus},
    struc::RegionLabel,
};
use bio::alignment::AlignmentOperation;
use itertools::Itertools;

pub type AlignOp = AlignmentOperation;

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

#[derive(Debug, PartialEq)]
pub enum Shape {
    Rect,
    HLine,
    VLine,
}

#[derive(Debug)]
pub struct PipeSeg {
    pub width: u32,
    pub color: Color,
    pub shape: Shape,
}

#[derive(Debug)]
pub struct Beta {
    pub value: f64,
    pub pos: usize,
}

#[derive(Debug)]
pub struct Pipe {
    pub segs: Vec<PipeSeg>,
    pub betas: Vec<Beta>,
    pub height: u32,
    pub outline: Vec<PipeSeg>,
    pub scale: Vec<(u32, Option<u32>)>,
}

pub struct Legend {
    pub labels: Vec<(String, Color)>,
    pub height: u32,
}

#[derive(Debug)]
pub struct DisplayParams {
    pub what_to_show: String,
    pub max_width: usize,
}

pub type PlotPanel = Vec<Pipe>;

pub struct PipePlot {
    pub panels: Vec<PlotPanel>,
    pub legend: Legend,
    // pub has_meth: bool,
}

pub fn encode_color(color: &Color) -> String {
    match color {
        Color::Purple => "#814ED1".to_string(),
        Color::Blue => "#1383C6".to_string(),
        Color::Orange => "#E16A2C".to_string(),
        Color::Teal => "#009CA2".to_string(),
        Color::Gray => "#BABABA".to_string(),
        Color::LightGray => "#D1D1D1".to_string(),
        Color::Black => "#000000".to_string(),
        Color::Pink => "#ED3981".to_string(),
        Color::Yellow => "#EFCD17".to_string(),
        Color::Green => "#009D4E".to_string(),
        Color::Red => "#E3371E".to_string(),
        Color::Khaki => "#F0E68C".to_string(),
        Color::PaleRed => "#FF4858".to_string(),
        Color::PaleBlue => "#46B2E8".to_string(),
        Color::Grad(value) => get_gradient(*value),
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

/*
pub fn get_base_pipe(locus: &Locus, allele: &Allele) -> Pipe {
    let mut segs = Vec::new();
    for label in &allele.all_labels {
        let color = get_color(&locus, AlignOp::Match, label);
        if let RegionLabel::Flank(start, end) = label {
            segs.push(PipeSeg {
                width: (*end - *start) as u32,
                color: color,
                shape: Shape::Rect,
            });
        } else if let RegionLabel::Seq(start, end) = label {
            segs.push(PipeSeg {
                width: (*end - *start) as u32,
                color: color,
                shape: Shape::Rect,
            });
        } else if let RegionLabel::Other(start, end) = label {
            segs.push(PipeSeg {
                width: (*end - *start) as u32,
                color: Color::Gray,
                shape: Shape::Rect,
            });
        } else if let RegionLabel::Tr(start, end, motif) = label {
            let tr_seq = &allele.seq[*start..*end];
            let mut index = 0;
            let mut last_matched_index = 0;
            while index as i32 <= tr_seq.len() as i32 - motif.len() as i32 {
                let observed_motif = &tr_seq[index..index + motif.len()];
                if observed_motif == motif {
                    if last_matched_index != index {
                        segs.push(PipeSeg {
                            width: (index - last_matched_index) as u32,
                            color: Color::Gray,
                            shape: Shape::Rect,
                        });
                    }
                    segs.push(PipeSeg {
                        width: motif.len() as u32,
                        color: color.clone(),
                        shape: Shape::Rect,
                    });
                    index += motif.len();
                    last_matched_index = index;
                } else {
                    index += 1;
                }
            }
            if last_matched_index != index {
                segs.push(PipeSeg {
                    width: (index - last_matched_index) as u32,
                    color: Color::Gray,
                    shape: Shape::Rect,
                });
            }
        } else {
            panic!();
        }
    }

    let height = 3;
    Pipe {
        segs,
        betas: Vec::new(),
        height,
        outline: Vec::new(),
        scale: Vec::new(),
    }
}
*/

pub fn get_allele_pipe(
    locus: &Locus,
    tick_spacing: usize,
    spans: &[RegionLabel],
    allele_labels: &Vec<BaseLabel>,
) -> Pipe {
    let mut span_and_label = Vec::new();
    let mut base_pos = 0;
    for label in allele_labels
        .iter()
        .filter(|l| **l != BaseLabel::MotifBound)
    {
        let span_index = get_label_index(spans, &AlignOp::Match, base_pos);
        span_and_label.push((span_index, label));

        base_pos += match label {
            BaseLabel::Match => 1,
            BaseLabel::Mismatch => 1,
            BaseLabel::MotifBound => 0,
            BaseLabel::NoMatch => 1,
            BaseLabel::Skip => 0,
        };
    }

    let mut segs = Vec::new();
    let groups = span_and_label.iter().group_by(|(a, b)| (a, b));
    for (key, group) in &groups {
        let length = group.count();
        let (span_index, allele_label) = key;

        let op = match allele_label {
            BaseLabel::Match => AlignOp::Match,
            BaseLabel::Mismatch => AlignOp::Subst,
            BaseLabel::MotifBound => AlignOp::Ins,
            BaseLabel::Skip => AlignOp::Ins,
            BaseLabel::NoMatch => AlignOp::Del,
        };

        let width = match op {
            AlignOp::Ins => 0,
            _ => length as u32,
        };

        let color = get_color(locus, op, &spans[*span_index]);

        let shape = match op {
            AlignOp::Del => Shape::HLine,
            AlignOp::Ins => Shape::VLine,
            _ => Shape::Rect,
        };

        segs.push(PipeSeg {
            width,
            color,
            shape,
        });
    }

    let outline = generate_outline(spans);
    let scale = generate_scale(tick_spacing, spans, allele_labels);

    Pipe {
        segs,
        betas: Vec::new(),
        height: 3,
        outline,
        scale,
    }
}

fn generate_scale(
    tick_spacing: usize,
    regions: &[RegionLabel],
    bases: &Vec<BaseLabel>,
) -> Vec<(u32, Option<u32>)> {
    let label_spacing = tick_spacing * 5;
    let mut bounds_by_tr = Vec::new();
    let mut allele_pos: u32 = 0;

    let mut motif_index = 0;
    for base in bases {
        if *base == BaseLabel::MotifBound {
            let region_index = get_region_id(regions, allele_pos);
            let is_first_motif = match &regions[region_index] {
                RegionLabel::Tr(start, _, _) => allele_pos == *start as u32,
                _ => continue,
            };

            if is_first_motif {
                motif_index = 0;
            } else {
                motif_index += 1;
            }

            if motif_index % tick_spacing != 0 {
                continue;
            }

            let label = if motif_index % label_spacing == 0 {
                Some(motif_index as u32)
            } else {
                None
            };

            bounds_by_tr.push((allele_pos, label));
        }

        allele_pos += match base {
            BaseLabel::Match => 1,
            BaseLabel::Mismatch => 1,
            BaseLabel::NoMatch => 1,
            _ => 0,
        };
    }

    bounds_by_tr
}

fn get_region_id(regions: &[RegionLabel], pos: u32) -> usize {
    for (index, region) in regions.iter().enumerate() {
        let (start, end) = match region {
            RegionLabel::Flank(start, end) => (start, end),
            RegionLabel::Other(start, end) => (start, end),
            RegionLabel::Seq(start, end) => (start, end),
            RegionLabel::Tr(start, end, _) => (start, end),
        };

        if *start as u32 <= pos && pos < *end as u32 {
            return index;
        }
    }

    panic!()
}

fn generate_outline(regions: &[RegionLabel]) -> Vec<PipeSeg> {
    regions
        .iter()
        .map(|r| {
            let width = match r {
                RegionLabel::Flank(start, end) => end - start,
                RegionLabel::Other(start, end) => end - start,
                RegionLabel::Seq(start, end) => end - start,
                RegionLabel::Tr(start, end, _) => end - start,
            } as u32;

            PipeSeg {
                width,
                color: Color::Black,
                shape: Shape::Rect,
            }
        })
        .collect_vec()
}

pub fn get_pipe(
    locus: &Locus,
    labels: &[RegionLabel],
    align: &AlignInfo,
    show_betas: bool,
) -> Pipe {
    assert!(align.align.xstart == 0);
    assert!(align.align.ystart == 0);

    let mut segs = Vec::new();
    let mut betas = Vec::new();

    let mut ops_and_locs = Vec::new();
    let mut allele_pos = 0;
    let mut seq_pos = 0;
    let mut cpg_index = 0;
    for op in &align.align.operations {
        let label_index = get_label_index(labels, op, allele_pos);
        ops_and_locs.push((*op, label_index));

        let is_cpg = seq_pos + 1 < align.seq.len() && &align.seq[seq_pos..seq_pos + 2] == "CG";

        if show_betas && *op != AlignOp::Del && align.meth.is_some() && is_cpg {
            if *op == AlignOp::Match || *op == AlignOp::Subst {
                let meth = align.meth.as_ref().unwrap();
                let value = meth[cpg_index] as f64 / 255.0;
                betas.push(Beta {
                    value,
                    pos: allele_pos,
                });
            }
            cpg_index += 1;
        }

        allele_pos += match *op {
            AlignOp::Match => 1,
            AlignOp::Subst => 1,
            AlignOp::Del => 1,
            AlignOp::Ins => 0,
            _ => panic!("Unhandled operation {:?}", *op),
        };

        seq_pos += match *op {
            AlignOp::Match => 1,
            AlignOp::Subst => 1,
            AlignOp::Del => 0,
            AlignOp::Ins => 1,
            _ => panic!("Unhandled operation {:?}", *op),
        };
    }

    assert_eq!(align.seq.len(), seq_pos);

    let groups = ops_and_locs.iter().group_by(|(a, l)| (a, l));
    for (key, group) in &groups {
        let group = group.collect_vec();
        segs.push(get_seg(locus, labels, key, group));
    }

    let height = 3;
    Pipe {
        segs,
        betas,
        height,
        outline: Vec::new(),
        scale: Vec::new(),
    }
}

// NOTE: ref_pos is not needed because we are working with global alignment
fn get_seg(
    locus: &Locus,
    labels: &[RegionLabel],
    op_and_seg_index: (&AlignOp, &usize),
    group: Vec<&(AlignOp, usize)>,
) -> PipeSeg {
    let (op, seg_index) = op_and_seg_index;
    let width = match op {
        AlignOp::Ins => 0,
        _ => group.len() as u32,
    };

    let color = get_color(locus, *op, &labels[*seg_index]);

    let shape = match op {
        AlignOp::Del => Shape::HLine,
        AlignOp::Ins => Shape::VLine,
        _ => Shape::Rect,
    };

    PipeSeg {
        width,
        color,
        shape,
    }
}

fn get_label_index(labels: &[RegionLabel], op: &AlignOp, pos: usize) -> usize {
    for (index, label) in labels.iter().enumerate() {
        let (start, end) = match label {
            RegionLabel::Flank(start, end) => (*start, *end),
            RegionLabel::Seq(start, end) => (*start, *end),
            RegionLabel::Tr(start, end, _) => (*start, *end),
            RegionLabel::Other(start, end) => (*start, *end),
        };

        if start <= pos && pos < end {
            return index;
        }

        if pos == end && *op == AlignOp::Ins {
            return index;
        }
    }

    println!("labels = {:?}, pos = {}", labels, pos);
    panic!();
}

pub fn get_color(locus: &Locus, op: AlignmentOperation, label: &RegionLabel) -> Color {
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

    if op == AlignOp::Subst {
        return Color::Gray;
    }

    if op == AlignOp::Ins {
        return Color::Black;
    }

    match label {
        RegionLabel::Flank(_, _) => Color::Teal,
        RegionLabel::Seq(_, _) => Color::LightGray,
        RegionLabel::Tr(_, _, motif) => {
            let index = locus.motifs.iter().position(|m| m == motif).unwrap();
            tr_colors[index % tr_colors.len()].clone()
        }
        RegionLabel::Other(_, _) => Color::LightGray,
    }
}
