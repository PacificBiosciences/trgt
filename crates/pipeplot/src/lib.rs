/*!
This crate provides functionality to generate the so-called "pipe plots"
consisting of stacked horizontal pipes (bars). Each pipe consists of segments of
specified width, shape, and color. Pipe plots can be annotated with scales,
labels, and legends. The crate supports rendering of pipe plots as SVG, PNG, and
PDF images.

Pipe plots are useful for representing pileups of sequenced reads.
*/

mod image;
mod pdf;
mod pipeplot;
mod png;
mod svg;

pub use image::generate as generate_image;
pub use pipeplot::{Band, Color, Legend, Pipe, PipePlot, Seg, Shape};
