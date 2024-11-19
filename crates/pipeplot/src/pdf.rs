use std::path::PathBuf;
use usvg::{TreeParsing, TreeTextToPath};

pub fn render(svg_path: &PathBuf, pdf_path: &PathBuf) -> Result<(), String> {
    let svg = std::fs::read_to_string(svg_path).map_err(|e| e.to_string())?;
    let options = usvg::Options::default();
    let mut tree = usvg::Tree::from_str(&svg, &options).expect("Error parsing SVG");
    let mut db = usvg::fontdb::Database::new();
    db.load_system_fonts();
    tree.convert_text(&db);
    let pdf = svg2pdf::convert_tree(&tree, svg2pdf::Options::default());
    std::fs::write(pdf_path, pdf).map_err(|e| e.to_string())?;
    Ok(())
}
