use std::path::PathBuf;
use usvg::{TreeParsing, TreeTextToPath};

pub fn render(svg_path: &PathBuf, png_path: &PathBuf) -> Result<(), String> {
    let svg = std::fs::read(svg_path).map_err(|e| e.to_string())?;
    let options = usvg::Options::default();
    let mut tree = usvg::Tree::from_data(&svg, &options).map_err(|e| e.to_string())?;
    let mut db = usvg::fontdb::Database::new();
    db.load_system_fonts();
    tree.convert_text(&db);
    let pixmap_size = tree.size.to_int_size();
    let mut pixmap = tiny_skia::Pixmap::new(pixmap_size.width(), pixmap_size.height())
        .ok_or("Unable to init image".to_string())?;
    let tree = resvg::Tree::from_usvg(&tree);
    tree.render(resvg::usvg::Transform::identity(), &mut pixmap.as_mut());
    pixmap.save_png(png_path).map_err(|e| e.to_string())?;
    Ok(())
}
