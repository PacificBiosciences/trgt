use crate::{pdf, png, svg, PipePlot};
use std::path::Path;

pub fn generate(plot: &PipePlot, path: &Path) -> Result<(), String> {
    if let Some(extension) = path.extension().and_then(|ext| ext.to_str()) {
        let svg_content = svg::generate_string(plot);
        match FileType::from_extension(extension) {
            Some(FileType::Svg) => svg::render_from_string(&svg_content, path),
            Some(FileType::Png) => png::render_from_string(&svg_content, path),
            Some(FileType::Pdf) => pdf::render_from_string(&svg_content, path),
            None => Err(format!("Unsupported file extension: {extension:?}")),
        }
    } else {
        Err(format!("Failed to get extension from path: {path:?}"))
    }
}

#[derive(Debug, PartialEq)]
enum FileType {
    Svg,
    Png,
    Pdf,
}

impl FileType {
    fn from_extension(extension: &str) -> Option<Self> {
        match extension.to_lowercase().as_str() {
            "svg" => Some(FileType::Svg),
            "png" => Some(FileType::Png),
            "pdf" => Some(FileType::Pdf),
            _ => None,
        }
    }
}
