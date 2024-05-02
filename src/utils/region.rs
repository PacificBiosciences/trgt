use crate::utils::Result;

#[derive(Debug, PartialEq)]
pub struct GenomicRegion {
    pub contig: String,
    pub start: u32,
    pub end: u32,
}

impl GenomicRegion {
    pub fn new(contig: impl Into<String>, start: u32, end: u32) -> Result<Self> {
        if start >= end {
            return Err(format!("Invalid region: start {} >= end {}", start, end));
        }

        Ok(Self {
            contig: contig.into(),
            start,
            end,
        })
    }

    pub fn from_string(encoding: &str) -> Result<Self> {
        let error_msg = || format!("Invalid region encoding: {}", encoding);
        let elements: Vec<&str> = encoding.split(&[':', '-']).collect();

        if elements.len() != 3 {
            return Err(error_msg());
        }

        let start: u32 = elements[1].parse().map_err(|_| error_msg())?;
        let end: u32 = elements[2].parse().map_err(|_| error_msg())?;

        Self::new(elements[0].to_string(), start, end)
    }

    pub fn intersect_position(&self, position: u32) -> bool {
        position >= self.start && position <= self.end
    }
}

#[cfg(test)]
mod tests {
    use super::GenomicRegion;
    #[test]
    fn init_region_from_valid_string_ok() {
        let region = GenomicRegion::from_string("chr1:100-200").unwrap();
        assert_eq!(region.contig, "chr1");
        assert_eq!(region.start, 100);
        assert_eq!(region.end, 200);
    }

    #[test]
    fn init_region_from_invalid_string_err() {
        assert_eq!(
            GenomicRegion::from_string("chr:1:100-200"),
            Err("Invalid region encoding: chr:1:100-200".to_string())
        );
    }

    #[test]
    fn init_region_from_invalid_start_err() {
        assert_eq!(
            GenomicRegion::from_string("chr:1:a-200"),
            Err("Invalid region encoding: chr:1:a-200".to_string())
        );
    }

    #[test]
    fn init_region_from_invalid_interval_err() {
        assert_eq!(
            GenomicRegion::from_string("chr1:200-100"),
            Err("Invalid region: start 200 >= end 100".to_string())
        );
    }

    #[test]
    fn init_region_from_invalid_interval_new() {
        assert_eq!(
            GenomicRegion::new("chr1", 200, 100),
            Err("Invalid region: start 200 >= end 100".to_string())
        );
    }
}
