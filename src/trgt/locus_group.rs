use crate::trgt::{
    genotype_reader::BamRead,
    locus::Locus,
    reads::{get_rq_tag, HiFiRead},
    workflows::Params,
};
use crate::utils::{GenomicRegion, Result};
use rand::{rngs::StdRng, Rng, SeedableRng};
use rust_htslib::bam::Record;
use std::sync::{Arc, Once};

static UNSORTED_WARNING: Once = Once::new();

/// A group of loci that are spatially close on the same contig.
#[derive(Debug)]
pub struct LocusGroup {
    loci: Vec<Locus>,
    pub region: GenomicRegion,
}

impl LocusGroup {
    /// Creates a new group with the first locus.
    pub fn new(locus: Locus) -> Self {
        let region = locus.region.clone();
        LocusGroup {
            loci: vec![locus],
            region,
        }
    }

    /// Checks if a locus can be added to the group.
    #[inline]
    pub fn can_add_locus(&self, locus: &Locus, max_span: u32) -> bool {
        if self.region.contig != locus.region.contig {
            return false;
        }

        // Check if the new locus starts after or at the last locus in the group.
        if locus.region.start < self.last_locus().region.start {
            UNSORTED_WARNING.call_once(|| {
                log::warn!(
                    "The input TR catalog appears to be unsorted! Consider sorting it with 'sort -k1,1 -k2,2n' for better performance. If this is not possible consider using more fetcher threads with the --fetcher-threads flag"
                );
            });
            return false;
        }

        // Check if the new locus would make the group span too large
        let new_end = self.region.end.max(locus.region.end);
        let new_span = new_end.saturating_sub(self.region.start);
        if new_span > max_span {
            return false;
        }
        true
    }

    /// Adds a locus to the group without checks, extending the group's region
    #[inline]
    pub fn add_locus_unchecked(&mut self, locus: Locus) {
        self.region.end = self.region.end.max(locus.region.end);
        self.loci.push(locus);
    }

    #[inline(always)]
    pub fn first_locus(&self) -> &Locus {
        self.loci.first().expect("non-empty group")
    }

    #[inline(always)]
    pub fn last_locus(&self) -> &Locus {
        self.loci.last().expect("non-empty group")
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.loci.len()
    }

    /// Checks whether the group contains no loci.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.loci.is_empty()
    }
}

#[inline]
fn locus_search_start(locus: &Locus, flank_len: u32) -> u32 {
    locus.region.start.saturating_sub(flank_len)
}

#[inline]
fn locus_search_end(locus: &Locus, flank_len: u32) -> u32 {
    locus.region.end.saturating_add(flank_len)
}

pub fn extract_reads_for_group<R: BamRead>(
    reader: &mut R,
    mut group: LocusGroup,
    params: &Arc<Params>,
) -> Result<Vec<Locus>> {
    let flank_len = params.search_flank_len as u32;
    let clip_radius = 2 * params.search_flank_len as i64;
    let group_region = &group.region;
    let reservoir_threshold = params.max_depth * 3;

    let extraction_region = (
        group_region.contig.as_str(),
        group_region.start.saturating_sub(flank_len),
        group_region.end.saturating_add(flank_len),
    );

    if let Err(msg) = reader.fetch(extraction_region) {
        return Err(format!(
            "BAM fetch error for group on {}: {}",
            group_region.contig, msg
        ));
    }

    let mut rng = StdRng::seed_from_u64(42);
    let mut n_reads_per_locus = vec![0; group.loci.len()];
    let mut locus_idx = 0;
    let mut record = Record::new();
    while let Some(result) = reader.read(&mut record) {
        match result {
            Ok(_) => {
                if record.is_supplementary() || record.is_secondary() || record.is_unmapped() {
                    continue;
                }
                if get_rq_tag(&record).unwrap_or(1.0) < params.min_read_qual {
                    continue;
                }

                let read_start = record.pos() as u32;
                let read_end = record.cigar().end_pos() as u32;

                // Find the first locus that could possibly overlap with the current read. Since both reads and loci are sorted by position, the locus pointer can be advanced
                while locus_idx < group.loci.len()
                    && locus_search_end(&group.loci[locus_idx], flank_len) < read_start
                {
                    locus_idx += 1;
                }

                let mut hifi_read_opt = None;
                for (i, locus) in group.loci.iter_mut().enumerate().skip(locus_idx) {
                    let locus_search_start = locus_search_start(locus, flank_len);
                    let locus_search_end = locus_search_end(locus, flank_len);
                    // The loci are sorted, if the locus start is past the read end, no subsequent loci can overlap
                    if locus_search_start > read_end {
                        break;
                    }
                    if read_end > locus_search_start && read_start < locus_search_end {
                        if hifi_read_opt.is_none() {
                            hifi_read_opt = Some(HiFiRead::from_hts_rec(&record));
                        }

                        if let Some(hifi_read) = &hifi_read_opt {
                            let region_to_clip = (
                                locus.region.start as i64 - clip_radius,
                                locus.region.end as i64 + clip_radius,
                            );

                            let current_count = n_reads_per_locus[i];
                            n_reads_per_locus[i] += 1;

                            if current_count < reservoir_threshold {
                                if let Some(clipped_read) = hifi_read.clip_to_region(region_to_clip)
                                {
                                    locus.reads.push(clipped_read);
                                }
                            } else {
                                let j = rng.random_range(0..=current_count);
                                if j < reservoir_threshold {
                                    if let Some(clipped_read) =
                                        hifi_read.clip_to_region(region_to_clip)
                                    {
                                        if j < locus.reads.len() {
                                            locus.reads[j] = clipped_read;
                                        } else {
                                            locus.reads.push(clipped_read);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            Err(_) => result.map_err(|e| e.to_string())?,
        }
    }

    for (i, locus) in group.loci.iter().enumerate() {
        let total_reads = n_reads_per_locus[i];
        if total_reads > reservoir_threshold {
            log::warn!(
                "{}: Reservoir sampled {} out of {} reads",
                locus.id,
                locus.reads.len(),
                total_reads
            );
        } else {
            log::debug!("{}: Collected {} reads", locus.id, locus.reads.len());
        }
    }

    Ok(group.loci)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::{Genotyper, Ploidy};

    fn create_test_locus(id: &str, contig: &str, start: u32, end: u32) -> Locus {
        Locus {
            id: id.to_string(),
            left_flank: b"AAAA".to_vec(),
            tr: b"AAAA".to_vec(),
            right_flank: b"AAAA".to_vec(),
            region: GenomicRegion {
                contig: contig.to_string(),
                start,
                end,
            },
            motifs: vec![b"A".to_vec()],
            struc: "(A)n".to_string(),
            ploidy: Ploidy::Two,
            genotyper: Genotyper::Size,
            reads: Vec::new(),
        }
    }

    #[test]
    fn test_new_group() {
        let locus = create_test_locus("test1", "chr1", 100, 200);
        let group = LocusGroup::new(locus);

        assert_eq!(group.len(), 1);
        assert!(!group.is_empty());
        assert_eq!(group.region.contig, "chr1");
        assert_eq!(group.region.start, 100);
        assert_eq!(group.region.end, 200);
        assert_eq!(group.first_locus().id, "test1");
        assert_eq!(group.last_locus().id, "test1");
    }

    #[test]
    fn test_can_add_locus_same_contig() {
        let locus1 = create_test_locus("test1", "chr1", 100, 200);
        let mut group = LocusGroup::new(locus1);
        let locus2 = create_test_locus("test2", "chr1", 250, 350);

        assert!(group.can_add_locus(&locus2, 1000));
        group.add_locus_unchecked(locus2);

        assert_eq!(group.len(), 2);
        assert_eq!(group.region.end, 350);
    }

    #[test]
    fn test_can_add_locus_different_contig() {
        let locus1 = create_test_locus("test1", "chr1", 100, 200);
        let group = LocusGroup::new(locus1);
        let locus2 = create_test_locus("test2", "chr2", 100, 200);

        assert!(!group.can_add_locus(&locus2, 1000));
    }

    #[test]
    fn test_can_add_locus_exceeds_max_span() {
        let locus1 = create_test_locus("test1", "chr1", 100, 200);
        let group = LocusGroup::new(locus1);
        let locus2 = create_test_locus("test2", "chr1", 600, 700);

        assert!(!group.can_add_locus(&locus2, 500)); // span would be 600
        assert!(group.can_add_locus(&locus2, 700)); // span would be 600, within limit
    }

    #[test]
    fn test_can_add_locus_unsorted_warning() {
        let locus1 = create_test_locus("test1", "chr1", 200, 300);
        let group = LocusGroup::new(locus1);
        let locus2 = create_test_locus("test2", "chr1", 100, 150); // starts before locus1

        assert!(!group.can_add_locus(&locus2, 1000));
    }

    #[test]
    fn test_add_locus_unchecked_extends_region() {
        let locus1 = create_test_locus("test1", "chr1", 100, 200);
        let mut group = LocusGroup::new(locus1);
        let locus2 = create_test_locus("test2", "chr1", 150, 300);

        assert!(group.can_add_locus(&locus2, 500));
        group.add_locus_unchecked(locus2);

        assert_eq!(group.len(), 2);
        assert_eq!(group.region.start, 100);
        assert_eq!(group.region.end, 300);
        assert_eq!(group.last_locus().id, "test2");
    }

    #[test]
    fn test_add_locus_unchecked_no_extension_needed() {
        let locus1 = create_test_locus("test1", "chr1", 100, 300);
        let mut group = LocusGroup::new(locus1);
        let locus2 = create_test_locus("test2", "chr1", 150, 250); // contained in locus1

        assert!(group.can_add_locus(&locus2, 500));
        group.add_locus_unchecked(locus2);

        assert_eq!(group.len(), 2);
        assert_eq!(group.region.start, 100);
        assert_eq!(group.region.end, 300);
    }

    #[test]
    fn test_multiple_loci_ordering() {
        let locus1 = create_test_locus("test1", "chr1", 100, 200);
        let mut group = LocusGroup::new(locus1);
        let locus2 = create_test_locus("test2", "chr1", 250, 350);
        let locus3 = create_test_locus("test3", "chr1", 400, 500);

        assert!(group.can_add_locus(&locus2, 500));
        group.add_locus_unchecked(locus2);
        assert!(group.can_add_locus(&locus3, 500));
        group.add_locus_unchecked(locus3);

        assert_eq!(group.len(), 3);
        assert_eq!(group.first_locus().id, "test1");
        assert_eq!(group.last_locus().id, "test3");
        assert_eq!(group.region.start, 100);
        assert_eq!(group.region.end, 500);
    }

    #[test]
    fn test_span_calculation_edge_cases() {
        let locus1 = create_test_locus("test1", "chr1", 0, 100);
        let group = LocusGroup::new(locus1);

        let locus2 = create_test_locus("test2", "chr1", 1000, 2000);
        let can_add = group.can_add_locus(&locus2, u32::MAX);
        assert!(can_add);

        let can_add_small_span = group.can_add_locus(&locus2, 500);
        assert!(!can_add_small_span);
    }

    #[test]
    fn test_exact_max_span_boundary() {
        let locus1 = create_test_locus("test1", "chr1", 100, 200);
        let group = LocusGroup::new(locus1);

        let locus2 = create_test_locus("test2", "chr1", 300, 400);

        assert!(group.can_add_locus(&locus2, 300));
        assert!(!group.can_add_locus(&locus2, 299));
    }

    #[test]
    fn test_is_empty_reflects_loci_storage() {
        let locus = create_test_locus("test1", "chr1", 100, 200);
        let mut group = LocusGroup::new(locus);

        group.loci.clear();

        assert!(group.is_empty());
    }

    #[test]
    fn test_locus_search_window_saturates_at_coordinate_bounds() {
        let locus = create_test_locus("test1", "chr1", 3, u32::MAX - 1);

        assert_eq!(locus_search_start(&locus, 10), 0);
        assert_eq!(locus_search_end(&locus, 10), u32::MAX);
    }
}
