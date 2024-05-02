#[derive(Debug, PartialEq, Clone)]
pub enum RegionLabel {
    Flank(usize, usize),      // Coordinates
    Tr(usize, usize, String), // Coordinates, Motif
    Seq(usize, usize),
    Other(usize, usize),
}
