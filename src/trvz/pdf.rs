use usvg::{TreeParsing, TreeTextToPath};

pub fn render(svg_path: &str, pdf_path: &str) {
    let svg = std::fs::read_to_string(svg_path).unwrap();
    let options = usvg::Options::default();
    let mut tree = usvg::Tree::from_str(&svg, &options).expect("Error parsing SVG");
    let mut db = usvg::fontdb::Database::new();
    db.load_system_fonts();
    tree.convert_text(&db);
    let pdf = svg2pdf::convert_tree(&tree, svg2pdf::Options::default());
    std::fs::write(pdf_path, pdf).unwrap();
}
