pub struct Read {
    pub seq: String,
    pub left_flank: usize,
    pub right_flank: usize,
    pub allele: i32,
    pub meth: Option<Vec<u8>>,
}
