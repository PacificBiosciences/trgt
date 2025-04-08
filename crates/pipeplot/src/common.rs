use usvg::Tree;

pub fn prepare_svg_tree(svg_data: &[u8]) -> Result<Tree, String> {
    let mut options = usvg::Options::default();
    let db = options.fontdb_mut();
    db.load_font_data(include_bytes!("../assets/fonts/RobotoMono-Bold.ttf").to_vec());
    db.load_system_fonts();
    let tree = usvg::Tree::from_data(svg_data, &options).map_err(|e| e.to_string())?;
    Ok(tree)
}
