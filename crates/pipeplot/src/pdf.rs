use crate::prepare_svg_tree;
use std::path::Path;

pub fn render_from_string(svg_content: &str, path: &Path) -> Result<(), String> {
    let tree = prepare_svg_tree(svg_content.as_bytes())?;
    let pdf = svg2pdf::to_pdf(
        &tree,
        svg2pdf::ConversionOptions::default(),
        svg2pdf::PageOptions::default(),
    )
    .map_err(|e| format!("Failed to convert SVG to PDF: {}", e))?;
    std::fs::write(path, pdf).map_err(|e| e.to_string())?;
    Ok(())
}
