use crate::pipeplot::{Color, FontConfig, Legend, Pipe, PipePlot, Shape};
use std::{fs, path::Path};

const DEFAULT_X_SCALE: f64 = 750.0;
const DEFAULT_Y_SCALE: f64 = 3.0;
const DEFAULT_PADDING: f64 = 12.0;

pub fn generate_string(plot: &PipePlot) -> String {
    let longest_pipe = get_longest_pipe(plot);
    let x_scale = DEFAULT_X_SCALE / longest_pipe as f64;
    let scale = (x_scale, DEFAULT_Y_SCALE);
    let mut generator = Generator::new(scale, DEFAULT_PADDING);
    generator.generate(plot);
    generator.buffer
}

pub fn render_from_string(svg_content: &str, path: &Path) -> Result<(), String> {
    fs::write(path, svg_content).map_err(|e| e.to_string())
}

fn get_longest_pipe(plot: &PipePlot) -> u32 {
    plot.pipes
        .iter()
        .map(|pipe| pipe.segs.iter().map(|seg| seg.width).sum())
        .max()
        .unwrap_or(0)
}

struct Generator {
    scale: (f64, f64),
    pad: f64,
    buffer: String,
}

impl Generator {
    fn new(scale: (f64, f64), pad: f64) -> Self {
        Self {
            scale,
            pad,
            buffer: String::with_capacity(10_000),
        }
    }

    pub fn generate(&mut self, pipe_plot: &PipePlot) {
        let (width, height) = self.get_dimensions(pipe_plot);
        self.start_svg(width, height);
        self.add_background();

        for pipe in &pipe_plot.pipes {
            self.plot_pipe(pipe, &pipe_plot.font);
            if pipe.outline {
                self.plot_outline(pipe);
            }
        }
        self.plot_legend(&pipe_plot.legend, &pipe_plot.font);
        self.end_svg();
    }

    fn add_line(&mut self, line: &str) {
        self.buffer.reserve(line.len() + 1);
        self.buffer.push_str(line);
        self.buffer.push('\n');
    }

    fn plot_legend(&mut self, legend: &Legend, font: &FontConfig) {
        let base_x = self.to_x(legend.xpos) + self.pad;
        let base_y = self.to_y(legend.ypos) + self.pad;
        let height = self.to_y(legend.height);
        let mut x = base_x;
        for (label, color) in &legend.labels {
            self.add_rect((x, base_y), (height, height), color, false);
            x += height + 2.0;
            let point = format!("x=\"{}\" y=\"{}\"", x, base_y + height - 1.0);
            let font_style = format!(
                r#"font-family="{}" font-weight="{}" font-size="{}""#,
                font.family, font.weight, font.size
            );
            let line = format!("<text {} {} >{}</text>", point, font_style, label);
            self.add_line(&line);
            x += 5.0 * (2 * label.len() as u32 + 1) as f64;
        }
    }

    fn plot_pipe(&mut self, pipe: &Pipe, font: &FontConfig) {
        let x = self.to_x(pipe.xpos) + self.pad;
        let y = self.to_y(pipe.ypos) + self.pad;
        let pipe_height = self.to_y(pipe.height);

        // Plot the main pipe
        let mut x_cur = x;

        for seg in &pipe.segs {
            let dims = (self.to_x(seg.width), pipe_height);
            match seg.shape {
                Shape::Rect => self.add_rect((x_cur, y), dims, &seg.color, true),
                Shape::HLine => self.add_hline((x_cur, y), dims, &seg.color, 1.5),
                Shape::Tick(label) => self.add_tick((x_cur, y), dims, &seg.color, label, font),
                Shape::None | Shape::VLine => {}
            }

            x_cur += self.to_x(seg.width);
        }

        // Plot insertions
        let mut x_cur = x;
        for seg in &pipe.segs {
            let dims = (self.to_x(seg.width), pipe_height);

            if seg.shape == Shape::VLine {
                self.add_vline((x_cur, y), dims, &seg.color);
            }

            x_cur += self.to_x(seg.width);
        }

        // Plot bands
        for band in &pipe.bands {
            let beta_x = x + self.to_x(band.pos);
            let dims = (self.to_x(1), pipe_height);
            self.add_rect((beta_x, y), dims, &band.color, false);
        }
    }

    fn plot_outline(&mut self, pipe: &Pipe) {
        let height = self.to_y(pipe.height);
        let width = self.to_x(pipe.segs.iter().map(|seg| seg.width).sum());

        let x = self.to_x(pipe.xpos) + self.pad;
        let y = self.to_y(pipe.ypos) + self.pad;

        let dimensions = format!("width=\"{width}\" height=\"{height}\"");
        let pos = format!("x=\"{x}\" y=\"{y}\"");
        let style = r##"stroke="#000000" stroke-width="1.5" fill="transparent""##;
        let line = format!("<rect {} {} {} />", dimensions, pos, style);
        self.add_line(&line);
    }

    fn add_rect(&mut self, pos: (f64, f64), dims: (f64, f64), color: &Color, add_highlight: bool) {
        let (x, y) = pos;
        let (w, h) = dims;

        let pos = format!("x=\"{}\" y=\"{}\"", x, y);
        let dim = format!("height=\"{}\" width=\"{}\"", h, w);

        let style = format!("fill=\"{}\" stroke=\"{}\" stroke-width=\"0\"", color, color);

        let rect = format!("<rect {} {} {} opacity=\"0.9\" />", pos, dim, style);
        self.add_line(&rect);

        if add_highlight {
            let pos = format!("x=\"{}\" y=\"{}\"", x, y + h * 0.18);
            let dim = format!("height=\"{}\" width=\"{}\"", h / 3.0, w);
            let highlight = format!("<rect {} {} fill=\"#F4EDF2\" opacity=\"0.25\" />", pos, dim);
            self.add_line(&highlight);
        }
    }

    fn add_hline(&mut self, pos: (f64, f64), dims: (f64, f64), color: &Color, stroke: f64) {
        let x1 = pos.0;
        let x2 = pos.0 + dims.0;
        let y1 = pos.1 + dims.1 / 2.0;
        let y2 = y1;

        let x1y1 = format!("x1=\"{}\" y1=\"{}\"", x1, y1);
        let x2y2 = format!("x2=\"{}\" y2=\"{}\"", x2, y2);

        let style = format!("stroke=\"{}\" stroke-width=\"{}\"", color, stroke);

        let line = format!("<line {} {} {} />", x1y1, x2y2, style);
        self.add_line(&line);
    }

    fn add_vline(&mut self, pos: (f64, f64), dims: (f64, f64), color: &Color) {
        let x1 = pos.0;
        let x2 = pos.0;
        let y1 = pos.1;
        let y2 = pos.1 + dims.1;

        let x1y1 = format!("x1=\"{}\" y1=\"{}\"", x1, y1);
        let x2y2 = format!("x2=\"{}\" y2=\"{}\"", x2, y2);

        let stroke_width = 2.0_f64.min(self.to_x(1));
        let style = format!("stroke=\"{}\" stroke-width=\"{}\"", color, stroke_width);

        let line = format!("<line {} {} {} />", x1y1, x2y2, style);
        self.add_line(&line);
    }

    fn add_tick(
        &mut self,
        pos: (f64, f64),
        dims: (f64, f64),
        color: &Color,
        label: Option<u32>,
        font: &FontConfig,
    ) {
        let x1 = pos.0;
        let x2 = pos.0;
        let y1 = pos.1;
        let y2 = pos.1 + dims.1;

        let x1y1 = format!("x1=\"{}\" y1=\"{}\"", x1, y1);
        let x2y2 = format!("x2=\"{}\" y2=\"{}\"", x2, y2);

        let stroke_width = 1.5;
        let style = format!("stroke=\"{}\" stroke-width=\"{}\"", color, stroke_width);

        let line = format!("<line {} {} {} />", x1y1, x2y2, style);
        self.add_line(&line);

        if let Some(label) = label {
            let point = format!("x=\"{}\" y=\"{}\"", x1, y1 - 2.0); // 2.0 is padding
            let font_style = format!(
                r#"font-family="{}" font-weight="{}" font-size="{}" text-anchor="middle""#,
                font.family, font.weight, font.size
            );
            let line = format!("<text {} {} >{}</text>", point, font_style, label);
            self.add_line(&line);
        }
    }

    fn get_dimensions(&self, pipe_plot: &PipePlot) -> (f64, f64) {
        let width = pipe_plot
            .pipes
            .iter()
            .map(|p| p.xpos + p.segs.iter().map(|s| s.width).sum::<u32>())
            .max()
            .unwrap_or(0);
        let height = pipe_plot.legend.ypos + pipe_plot.legend.height;

        let xdim = self.to_x(width) + 2.0 * self.pad;
        let ydim = self.to_y(height) + 2.0 * self.pad;
        (xdim, ydim)
    }

    fn start_svg(&mut self, width: f64, height: f64) {
        self.add_line(r#"<?xml version="1.0"?>"#);
        let line = format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{}" height="{}">"#,
            width, height
        );
        self.add_line(&line);
    }

    fn end_svg(&mut self) {
        self.add_line("</svg>");
    }

    fn add_background(&mut self) {
        self.add_line(r#"<rect width="100%" height="100%" fill="white"/>"#);
    }

    fn to_x(&self, x: u32) -> f64 {
        x as f64 * self.scale.0
    }

    fn to_y(&self, y: u32) -> f64 {
        y as f64 * self.scale.1
    }
}
