use crate::utils::Result;
use std::fmt;

#[derive(Debug, PartialEq, Clone)]
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
        let elements: Vec<&str> = encoding.split(&[':', '-']).collect();
        if elements.len() != 3 {
            return Err(format!(
                "Invalid region encoding '{}': expected 'chrom:start-end'",
                encoding
            ));
        }
        Self::from_bed_fields(&elements)
    }

    pub fn from_bed_fields(fields: &[&str]) -> Result<Self> {
        if fields.len() < 3 {
            return Err(format!(
                "Invalid region encoding: {:?}. BED fields must contain at least chrom, start, and end columns", fields
            ));
        }

        let start: u32 = fields[1].parse().map_err(|_| {
            format!(
                "Invalid region encoding: {:?}. Invalid start position '{}': must be a non-negative integer",
                fields, fields[1]
            )
        })?;

        let end: u32 = fields[2].parse().map_err(|_| {
            format!(
                "Invalid region encoding: {:?}. Invalid end position '{}': must be a non-negative integer",
                fields, fields[2]
            )
        })?;

        Self::new(fields[0].to_string(), start, end)
    }

    pub fn intersect_position(&self, position: u32) -> bool {
        position >= self.start && position <= self.end
    }

    pub fn as_parts(&self) -> (&str, u32, u32) {
        (self.contig.as_str(), self.start, self.end)
    }
}

impl fmt::Display for GenomicRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.start, self.end)
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
    fn display_formats_as_region_label() {
        let region = GenomicRegion::new("chr1", 50, 200).unwrap();

        assert_eq!(region.to_string(), "chr1:50-200");
    }

    #[test]
    fn exposes_region_parts() {
        let region = GenomicRegion::new("chr1", 50, 200).unwrap();

        assert_eq!(region.as_parts(), ("chr1", 50, 200));
    }

    #[test]
    fn init_region_from_invalid_string_err() {
        assert_eq!(
            GenomicRegion::from_string("chr:1:100-200"),
            Err("Invalid region encoding 'chr:1:100-200': expected 'chrom:start-end'".to_string())
        );
    }

    #[test]
    fn init_region_from_invalid_start_err() {
        assert_eq!(
            GenomicRegion::from_string("chr:1:a-200"),
            Err("Invalid region encoding 'chr:1:a-200': expected 'chrom:start-end'".to_string())
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

    #[test]
    fn init_region_from_bed_fields_valid() {
        let fields = vec!["chr1", "100", "200"];
        let region = GenomicRegion::from_bed_fields(&fields).unwrap();
        assert_eq!(region.contig, "chr1");
        assert_eq!(region.start, 100);
        assert_eq!(region.end, 200);
    }

    #[test]
    fn init_region_from_bed_fields_padding_valid() {
        let fields = vec!["chr1", "100", "200"];
        let region = GenomicRegion::from_bed_fields(&fields).unwrap();
        assert_eq!(region.contig, "chr1");
        assert_eq!(region.start, 100);
        assert_eq!(region.end, 200);
    }

    #[test]
    fn init_region_from_bed_fields_invalid_start() {
        let fields = vec!["chr1", "wrong", "200"];
        assert!(GenomicRegion::from_bed_fields(&fields).is_err());
    }

    #[test]
    fn init_region_from_bed_fields_too_few_fields() {
        let fields = vec!["chr1", "100"];
        assert!(GenomicRegion::from_bed_fields(&fields).is_err());
    }
}
