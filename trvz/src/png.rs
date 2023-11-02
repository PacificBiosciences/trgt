pub fn render(svg_path: &str, png_path: &str) {
    // Source: https://github.com/RazrFalcon/resvg/blob/master/examples/minimal.rs
    let mut opt = usvg::Options::default();
    // opt.resources_dir = std::fs::canonicalize(svg_path)
    //     .ok()
    //     .and_then(|p| p.parent().map(|p| p.to_path_buf()));
    opt.fontdb.load_system_fonts();

    let svg_data = std::fs::read(svg_path).unwrap();
    let rtree = usvg::Tree::from_data(&svg_data, &opt.to_ref()).unwrap();

    let pixmap_size = rtree.svg_node().size.to_screen_size();
    let mut pixmap = tiny_skia::Pixmap::new(pixmap_size.width(), pixmap_size.height()).unwrap();
    resvg::render(
        &rtree,
        usvg::FitTo::Original,
        tiny_skia::Transform::default(),
        pixmap.as_mut(),
    )
    .unwrap();
    pixmap.save_png(png_path).unwrap();
}
