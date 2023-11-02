#[derive(Debug, PartialEq)]
pub struct GenomicRegion {
    pub contig: String,
    pub start: u32,
    pub end: u32,
}

impl GenomicRegion {
    pub fn new(encoding: &str) -> Result<Self, String> {
        let error_msg = || format!("Invalid region encoding: {}", encoding);
        let elements: Vec<&str> = encoding.split(&[':', '-']).collect();

        if elements.len() != 3 {
            return Err(error_msg());
        }

        let start: u32 = elements[1].parse().map_err(|_| error_msg())?;

        let end: u32 = elements[2].parse().map_err(|_| error_msg())?;

        if start >= end {
            return Err(format!("Invalid region: start {} >= end {}", start, end));
        }

        Ok(Self {
            contig: elements[0].to_string(),
            start,
            end,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::GenomicRegion;
    #[test]
    fn init_region_from_valid_string_ok() {
        let region = GenomicRegion::new("chr1:100-200").unwrap();
        assert_eq!(region.contig, "chr1");
        assert_eq!(region.start, 100);
        assert_eq!(region.end, 200);
    }

    #[test]
    fn init_region_from_invalid_string_err() {
        assert_eq!(
            GenomicRegion::new("chr:1:100-200"),
            Err("Invalid region encoding: chr:1:100-200".to_string())
        );
    }
}
