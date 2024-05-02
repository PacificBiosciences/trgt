use crate::trvz::pipe_plot::{encode_color, Color, Legend, Pipe, PipePlot, PlotPanel, Shape};
use std::{fs::File, io::Write};

pub fn generate(genotype_plot: &PipePlot, path: &str) {
    let start = (10.0, 14.0);

    let longest_pipe = get_longest_pipe(genotype_plot);
    let x_scale = 750.0 / longest_pipe as f64;

    let scale = (x_scale, 3.0);
    let base_spacing = 3.0;
    let read_spacing = 2.0;
    let allele_spacing = 20.0;
    let legend_spacing = 4.0;
    let file = File::create(path).unwrap();
    let mut generator = Generator {
        start,
        scale,
        base_spacing,
        read_spacing,
        allele_spacing,
        legend_spacing,
        file,
    };
    generator.generate(genotype_plot);
}

fn get_longest_pipe(plot: &PipePlot) -> u32 {
    let mut widths = Vec::new();
    for panel in &plot.panels {
        for pipe in panel {
            let mut width = 0;
            for seg in &pipe.segs {
                width += seg.width;
            }
            widths.push(width);
        }
    }
    if widths.is_empty() {
        0
    } else {
        *widths.iter().max().unwrap()
    }
}

struct Generator {
    start: (f64, f64),
    scale: (f64, f64),
    base_spacing: f64,
    read_spacing: f64,
    allele_spacing: f64,
    legend_spacing: f64,
    file: File,
}

impl Generator {
    pub fn generate(&mut self, genotype_plot: &PipePlot) {
        let (width, height) = self.get_dimensions(genotype_plot);
        self.start_svg(width, height);
        self.add_background();

        let base_x = self.start.0;
        let mut base_y = self.start.1;
        for (index, allele_plot) in genotype_plot.panels.iter().enumerate() {
            self.plot_allele(base_x, &mut base_y, allele_plot);
            if index + 1 != genotype_plot.panels.len() {
                base_y += self.allele_spacing;
            }
        }

        base_y += self.legend_spacing;
        self.plot_legend(base_x, base_y, &genotype_plot.legend);

        self.end_svg();
    }

    fn plot_allele(&mut self, base_x: f64, base_y: &mut f64, allele_plot: &PlotPanel) {
        for pipe in allele_plot {
            self.plot_pipe(pipe, base_x, *base_y);
            if !pipe.outline.is_empty() {
                self.plot_outline(pipe, base_x, *base_y);
                self.plot_scale(pipe, base_x, *base_y);
                *base_y += self.to_y(pipe.height) + self.base_spacing;
            } else {
                *base_y += self.to_y(pipe.height) + self.read_spacing;
            }
        }
    }

    fn plot_legend(&mut self, base_x: f64, base_y: f64, legend: &Legend) {
        let height = self.to_y(legend.height);
        let mut x = base_x;
        for (label, color) in &legend.labels {
            self.add_rect((x, base_y), (height, height), color, false);
            x += height + 2.0; //self.to_y(1);
            let point = format!("x=\"{}\" y=\"{}\"", x, base_y + height - 1.0);
            let height = r#"font-size="14px""#;
            let style = r#"font-family="monospace" font-weight="bold""#;
            let line = format!("<text {} {} {} >{}</text>", point, style, height, label);
            writeln!(self.file, "{}", line).unwrap();
            x += 5.0 * (2 * label.len() as u32 + 1) as f64;
        }
    }

    fn plot_pipe(&mut self, pipe: &Pipe, x: f64, y: f64) {
        let pipe_height = self.to_y(pipe.height);

        // Plot the main pipe
        let mut x_cur = x;

        self.add_left_cap((x_cur, y), pipe_height);
        for seg in &pipe.segs {
            let dims = (self.to_x(seg.width), pipe_height);
            if seg.shape == Shape::Rect {
                self.add_rect((x_cur, y), dims, &seg.color, true);
            }

            if seg.shape == Shape::HLine {
                self.add_hline((x_cur, y), dims, &seg.color, 1.5);
            }

            x_cur += self.to_x(seg.width);
        }
        self.add_right_cap((x_cur, y), pipe_height);

        // Plot insertions
        let mut x_cur = x;
        for seg in &pipe.segs {
            let dims = (self.to_x(seg.width), pipe_height);

            if seg.shape == Shape::VLine {
                self.add_vline((x_cur, y), dims, &seg.color);
            }

            x_cur += self.to_x(seg.width);
        }

        // Plot betas
        for beta in &pipe.betas {
            let beta_x = x + self.to_x(beta.pos as u32);
            let dims = (self.to_x(1), pipe_height);
            self.add_rect((beta_x, y), dims, &Color::Grad(beta.value), false);
        }
    }

    fn plot_outline(&mut self, pipe: &Pipe, x: f64, y: f64) {
        let height = self.to_y(pipe.height);
        let mut width = 0;
        for seg in &pipe.outline {
            width += seg.width;
        }
        let width = self.to_x(width);

        // Draw segment boundaries
        let mut bar_x = x + self.to_x(pipe.outline.first().unwrap().width);
        for segment in pipe.outline.iter().skip(1) {
            let point1 = format!("x1=\"{}\" y1=\"{}\"", bar_x, y);
            let point2 = format!("x2=\"{}\" y2=\"{}\"", bar_x, y + height);
            let style = r##"stroke="#000000" stroke-width="1.5""##;
            let line = format!("<line {} {} {} />", point1, point2, style);
            writeln!(self.file, "{}", line).unwrap();
            bar_x += self.to_x(segment.width);
        }

        // Draw top horizontal line
        let point1 = format!("x1=\"{}\" y1=\"{}\"", x, y);
        let point2 = format!("x2=\"{}\" y2=\"{}\"", x + width, y);
        let style = r##"stroke="#000000" stroke-width="1.5""##;
        let line = format!("<line {} {} {} />", point1, point2, style);
        writeln!(self.file, "{}", line).unwrap();

        // Draw bottom horizontal line
        let point1 = format!("x1=\"{}\" y1=\"{}\"", x, y + height);
        let point2 = format!("x2=\"{}\" y2=\"{}\"", x + width, y + height);
        let style = r##"stroke="#000000" stroke-width="1.5""##;
        let line = format!("<line {} {} {} />", point1, point2, style);
        writeln!(self.file, "{}", line).unwrap();

        // Draw left cap
        let point1 = format!("{} {}", x, y);
        let point2 = format!("{} {}", x - height / 2.0, y);
        let point3 = format!("{} {}", x - height / 2.0, y + height);
        let point4 = format!("{} {}", x, y + height);

        let path = format!("d=\"M {} C {}, {}, {}\"", point1, point2, point3, point4);
        let style = r##"stroke="#000000" stroke-width="1.5" fill="none" opacity="0.9""##;

        let line = format!("<path {} {} />", path, style);
        writeln!(self.file, "{}", line).unwrap();

        // Draw right cap
        let point1 = format!("{} {}", x + width, y);
        let point2 = format!("{} {}", x + width + height / 2.0, y);
        let point3 = format!("{} {}", x + width + height / 2.0, y + height);
        let point4 = format!("{} {}", x + width, y + height);

        let path = format!("d=\"M {} C {}, {}, {}\"", point1, point2, point3, point4);
        let style = r##"stroke="#000000" stroke-width="1.5" fill="none" opacity="0.9""##;

        let line = format!("<path {} {} />", path, style);
        writeln!(self.file, "{}", line).unwrap();
    }

    fn plot_scale(&mut self, pipe: &Pipe, x: f64, y: f64) {
        let height = self.to_y(pipe.height) / 2.0;

        // Draw segment boundaries
        for (x_pos, label) in pipe.scale.iter() {
            let tick_x = x + self.to_x(*x_pos);
            let point1 = format!("x1=\"{}\" y1=\"{}\"", tick_x, y - height / 2.0);
            let point2 = format!("x2=\"{}\" y2=\"{}\"", tick_x, y + height / 2.0);
            let style = r##"stroke="#000000" stroke-width="1.5""##;
            let line = format!("<line {} {} {} />", point1, point2, style);
            writeln!(self.file, "{}", line).unwrap();

            if let Some(label) = *label {
                let point = format!("x=\"{}\" y=\"{}\"", tick_x, y - height / 2.0 - 2.0); // 2.0 is padding
                let height = r#"font-size="10px""#;
                let style = r#"font-family="monospace" font-weight="bold" text-anchor="middle""#;
                let line = format!("<text {} {} {} >{}</text>", point, style, height, label);
                writeln!(self.file, "{}", line).unwrap();
            }
        }
    }

    fn add_rect(&mut self, pos: (f64, f64), dims: (f64, f64), color: &Color, add_highlight: bool) {
        let (x, y) = pos;
        let (w, h) = dims;

        let pos = format!("x=\"{}\" y=\"{}\"", x, y);
        let dim = format!("height=\"{}\" width=\"{}\"", h, w);

        let color = encode_color(color);
        let style = format!("fill=\"{}\" stroke=\"{}\" stroke-width=\"0\"", color, color);

        let rect = format!("<rect {} {} {} opacity=\"0.9\" />", pos, dim, style);
        writeln!(self.file, "{}", rect).unwrap();

        if add_highlight {
            let pos = format!("x=\"{}\" y=\"{}\"", x, y + h * 0.18);
            let dim = format!("height=\"{}\" width=\"{}\"", h / 3.0, w);
            let highlight = format!("<rect {} {} fill=\"#F4EDF2\" opacity=\"0.25\" />", pos, dim);
            writeln!(self.file, "{}", highlight).unwrap();
        }
    }

    fn add_hline(&mut self, pos: (f64, f64), dims: (f64, f64), color: &Color, stroke: f64) {
        let x1 = pos.0;
        let x2 = pos.0 + dims.0;
        let y1 = pos.1 + dims.1 / 2.0;
        let y2 = y1;

        let x1y1 = format!("x1=\"{}\" y1=\"{}\"", x1, y1);
        let x2y2 = format!("x2=\"{}\" y2=\"{}\"", x2, y2);

        let style = format!(
            "stroke=\"{}\" stroke-width=\"{}\"",
            encode_color(color),
            stroke
        );

        let line = format!("<line {} {} {} />", x1y1, x2y2, style);
        writeln!(self.file, "{}", line).unwrap();
    }

    fn add_vline(&mut self, pos: (f64, f64), dims: (f64, f64), color: &Color) {
        let x1 = pos.0;
        let x2 = pos.0;
        let y1 = pos.1;
        let y2 = pos.1 + dims.1;

        let x1y1 = format!("x1=\"{}\" y1=\"{}\"", x1, y1);
        let x2y2 = format!("x2=\"{}\" y2=\"{}\"", x2, y2);

        let stroke_width = 2.0_f64.min(self.to_x(1));
        let style = format!(
            "stroke=\"{}\" stroke-width=\"{}\"",
            encode_color(color),
            stroke_width
        );

        let line = format!("<line {} {} {} />", x1y1, x2y2, style);
        writeln!(self.file, "{}", line).unwrap();
    }

    fn add_left_cap(&mut self, pos: (f64, f64), height: f64) {
        let (x, y) = pos;
        let point1 = format!("{} {}", x, y);
        let point2 = format!("{} {}", x - height / 2.0, y);
        let point3 = format!("{} {}", x - height / 2.0, y + height);
        let point4 = format!("{} {}", x, y + height);

        let path = format!("d=\"M {} C {}, {}, {}\"", point1, point2, point3, point4);
        let fill = r##"fill="#009CA2""##;
        let opacity = r#"opacity="0.9""#;

        let line = format!("<path {} {} {} />", path, fill, opacity);
        writeln!(self.file, "{}", line).unwrap();
    }

    fn add_right_cap(&mut self, pos: (f64, f64), height: f64) {
        let (x, y) = pos;
        let point1 = format!("{} {}", x, y);
        let point2 = format!("{} {}", x + height / 2.0, y);
        let point3 = format!("{} {}", x + height / 2.0, y + height);
        let point4 = format!("{} {}", x, y + height);

        let path = format!("d=\"M {} C {}, {}, {}\"", point1, point2, point3, point4);
        let fill = "fill=\"#009CA2\"";

        let opacity = r#"opacity="0.9""#;

        let line = format!("<path {} {} {} />", path, fill, opacity);
        writeln!(self.file, "{}", line).unwrap();
    }

    fn get_dimensions(&self, genotype_plot: &PipePlot) -> (f64, f64) {
        let mut widths = Vec::new();
        let mut height = self.start.1;

        for allele_plot in &genotype_plot.panels {
            for (index, pipe) in allele_plot.iter().enumerate() {
                if index == 0 {
                    height += self.to_y(pipe.height) + self.base_spacing;
                } else {
                    height += self.to_y(pipe.height) + self.read_spacing;
                }

                let mut pipe_width = 2.0 * self.start.0;
                for seg in &pipe.segs {
                    pipe_width += self.to_x(seg.width);
                }
                widths.push(pipe_width);
            }
        }

        height += self.allele_spacing * (genotype_plot.panels.len() as f64 - 1.0);
        height +=
            self.legend_spacing + self.to_y(genotype_plot.legend.height) + self.legend_spacing;

        let width = *widths
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        (width, height)
    }

    fn start_svg(&mut self, width: f64, height: f64) {
        writeln!(self.file, r#"<?xml version="1.0"?>"#).unwrap();
        let line = r#"<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" "#;
        write!(self.file, "{}", line).unwrap();
        writeln!(self.file, "width=\"{}\" height=\"{}\">", width, height).unwrap();
    }

    fn end_svg(&mut self) {
        writeln!(self.file, "</svg>").unwrap();
    }

    fn add_background(&mut self) {
        writeln!(
            self.file,
            r#"<rect width="100%" height="100%" fill="white"/>"#
        )
        .unwrap();
    }

    fn to_x(&self, x: u32) -> f64 {
        x as f64 * self.scale.0
    }

    fn to_y(&self, y: u32) -> f64 {
        y as f64 * self.scale.1
    }
}
