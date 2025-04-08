use crate::prepare_svg_tree;
use std::path::Path;

pub fn render_from_string(svg_content: &str, path: &Path) -> Result<(), String> {
    let tree = prepare_svg_tree(svg_content.as_bytes())?;
    let pixmap_size = tree.size().to_int_size();
    let mut pixmap = resvg::tiny_skia::Pixmap::new(pixmap_size.width(), pixmap_size.height())
        .ok_or("Unable to init image".to_string())?;
    resvg::render(
        &tree,
        resvg::usvg::Transform::identity(),
        &mut pixmap.as_mut(),
    );
    pixmap.save_png(path).map_err(|e| e.to_string())?;
    Ok(())
}
