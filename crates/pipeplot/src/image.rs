use crate::{pdf, png, svg, PipePlot};
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;

pub fn generate(plot: &PipePlot, path: &PathBuf) -> Result<(), String> {
    let extension = Path::new(path)
        .extension()
        .and_then(|ext| ext.to_str())
        .ok_or_else(|| format!("Failed to get extension from path: {path:?}"))?;
    let file_type = FileType::from_extension(extension)
        .ok_or_else(|| format!("Unsupported file extension: {extension:?}"))?;

    let temp_dir = PathBuf::from(&path);
    let temp_dir = temp_dir
        .parent()
        .ok_or_else(|| format!("Invalid path: {path:?}"))?;
    let svg_temp_path = NamedTempFile::new_in(temp_dir)
        .map_err(|e| format!("Failed to create temporary file: {e}"))?
        .into_temp_path();

    svg::generate(plot, &svg_temp_path.to_path_buf());
    match file_type {
        FileType::Svg => svg_temp_path.persist(path).map_err(|e| e.to_string())?,
        FileType::Png => png::render(&svg_temp_path.to_path_buf(), path)?,
        FileType::Pdf => pdf::render(&svg_temp_path.to_path_buf(), path)?,
    }
    Ok(())
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
