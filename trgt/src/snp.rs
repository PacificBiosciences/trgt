use crate::reads::cigar::{Cigar, CigarOp};
use crate::reads::HiFiRead;
use crate::utils::GenomicRegion;

impl GenomicRegion {
    pub fn intersect_position(&self, position: u32) -> bool {
        position >= self.start && position <= self.end
    }
}

// Collect all positions in a read that are mismatches relative to the starting point of the alignment
// extract_snps_offset should be used instead
// 0-based indexing
#[allow(dead_code)]
pub fn extract_snps(cigar: &Cigar, region: &GenomicRegion) -> Vec<u32> {
    let mut mismatches = Vec::new();
    let mut start_ref = cigar.ref_pos as u32;
    for op in cigar.ops.iter() {
        match op {
            // Only mismatches not part of the region should be tracked
            CigarOp::Diff(n) if !region.intersect_position(start_ref) => {
                mismatches.extend(start_ref..start_ref + *n);
                start_ref += *n;
            }
            // Operations that consume the reference (including mismatches within the region)
            CigarOp::Match(n)
            | CigarOp::Diff(n)
            | CigarOp::Equal(n)
            | CigarOp::Del(n)
            | CigarOp::RefSkip(n) => {
                start_ref += *n;
            }
            _ => {}
        }
    }
    mismatches
}

// Collect all positions in a read that are mismatches relative to the starting point of the region
#[allow(dead_code)]
pub fn extract_snps_offset(cigar: &Cigar, region: &GenomicRegion) -> Vec<i32> {
    let mut mismatches: Vec<i32> = Vec::new();
    let mut start_ref = cigar.ref_pos as u32;
    for op in cigar.ops.iter() {
        match op {
            // Only mismatches not part of the region should be tracked
            CigarOp::Diff(n) if !region.intersect_position(start_ref) => {
                let diff = if start_ref < region.start {
                    (start_ref as i32) - (region.start as i32)
                } else {
                    (start_ref as i32) - (region.end as i32)
                };
                mismatches.extend((0..*n).map(|i| diff + i as i32));
                start_ref += *n;
            }
            // Operations that consume the reference (including mismatches within the region)
            CigarOp::Match(n)
            | CigarOp::Diff(n)
            | CigarOp::Equal(n)
            | CigarOp::Del(n)
            | CigarOp::RefSkip(n) => {
                start_ref += *n;
            }
            _ => {}
        }
    }
    mismatches
}

#[allow(dead_code)]
pub fn analyze_snps(reads: &mut [HiFiRead], region: &GenomicRegion) {
    for read in reads.iter_mut().filter(|r| r.cigar.is_some()) {
        let mismatches = extract_snps_offset(read.cigar.as_ref().unwrap(), region);
        read.mismatch_offsets = Some(mismatches);
    }
}

#[allow(dead_code)]
fn get_cigar_string(cigar: &Vec<CigarOp>) -> String {
    let mut result = String::new();
    for op in cigar {
        result.push_str(&op.len().to_string());
        result.push(op.char());
    }
    result
}

#[allow(dead_code)]
fn parse_cigar_string(cigar_string: &str) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut num_str = String::new();
    for c in cigar_string.chars() {
        if c.is_ascii_digit() {
            num_str.push(c);
        } else {
            let num = if num_str.is_empty() {
                1
            } else {
                num_str.parse().unwrap()
            };
            num_str.clear();
            let op = match c {
                'M' => CigarOp::Match(num),
                'I' => CigarOp::Ins(num),
                'D' => CigarOp::Del(num),
                'N' => CigarOp::RefSkip(num),
                'S' => CigarOp::SoftClip(num),
                'H' => CigarOp::HardClip(num),
                'P' => CigarOp::Pad(num),
                '=' => CigarOp::Equal(num),
                'X' => CigarOp::Diff(num),
                _ => panic!("Invalid cigar operation: {}", c),
            };
            ops.push(op);
        }
    }
    ops
}

#[cfg(test)]
mod tests {
    use super::*;

    // Corresponds to: m64015_190920_185703/70518060/ccs in the example data of TRGT
    const CIGAR_STRING: &str = "178=1D6=1D51=1I59=1I42=1I22=1X1=2D7=2X22=1D82=1X66=1D8=3D261=1D20=1D8=2D60=1I116=1I9=2I5=1I70=1D4=1D3=1D20=1I25=1I30=1X18=7I205=1I85=1D230=1I25=1I191=1I66=1I116=1D25=1I18=1I135=1I88=1I15=1I33=1I139=1D11=1I81=1I100=1I131=1I37=1I20=1X94=1I9=1D152=1I8=1D33=1D115=1I198=1I8=1I37=1I14=1I60=1X88=1D203=1D16=1I24=1I153=1D97=1I108=1I9=2I80=1D23=1I12=1I171=1X19=1D75=1I31=1X3=1I12=1I15=1X47=1I104=2I106=1I4=1I32=1D20=1I3=1X61=1I15=1X30=1D118=1D86=1X28=1I40=1X34=1I25=1I7=1I7=1I36=1D100=1I487=1I340=1I198=1I5=1I95=1I203=1I59=1I3=1D16=1I30=1I26=1I121=1X114=1I150=1I215=2I209=1I27=1I7=1D14=1I410=1I87=1I404=1I19=1I3=1I99=1D525=27D178=1I177=1I326=1X32=1I13=1I54=1D33=1I87=1X39=1D90=3285S";
    #[test]
    fn test_cigar_conversion() {
        let cigar = Cigar {
            ref_pos: 372,
            ops: vec![
                CigarOp::Equal(178),
                CigarOp::Del(1),
                CigarOp::Equal(6),
                CigarOp::Del(1),
                CigarOp::Equal(51),
                CigarOp::Ins(1),
                CigarOp::Equal(59),
                CigarOp::Ins(1),
                CigarOp::Equal(42),
                CigarOp::Ins(1),
                CigarOp::Equal(22),
                CigarOp::Diff(1),
                CigarOp::Equal(1),
                CigarOp::Del(2),
                CigarOp::Equal(7),
                CigarOp::Diff(2),
                CigarOp::Equal(22),
                CigarOp::Del(1),
                CigarOp::Equal(82),
                CigarOp::Diff(1),
                CigarOp::Equal(66),
                CigarOp::Del(1),
                CigarOp::Equal(8),
                CigarOp::Del(3),
                CigarOp::Equal(261),
                CigarOp::Del(1),
                CigarOp::Equal(20),
                CigarOp::Del(1),
                CigarOp::Equal(8),
                CigarOp::Del(2),
                CigarOp::Equal(60),
                CigarOp::Ins(1),
                CigarOp::Equal(116),
                CigarOp::Ins(1),
                CigarOp::Equal(9),
                CigarOp::Ins(2),
                CigarOp::Equal(5),
                CigarOp::Ins(1),
                CigarOp::Equal(70),
                CigarOp::Del(1),
                CigarOp::Equal(4),
                CigarOp::Del(1),
                CigarOp::Equal(3),
                CigarOp::Del(1),
                CigarOp::Equal(20),
                CigarOp::Ins(1),
                CigarOp::Equal(25),
                CigarOp::Ins(1),
                CigarOp::Equal(30),
                CigarOp::Diff(1),
                CigarOp::Equal(18),
                CigarOp::Ins(7),
                CigarOp::Equal(205),
                CigarOp::Ins(1),
                CigarOp::Equal(85),
                CigarOp::Del(1),
                CigarOp::Equal(230),
                CigarOp::Ins(1),
                CigarOp::Equal(25),
                CigarOp::Ins(1),
                CigarOp::Equal(191),
                CigarOp::Ins(1),
                CigarOp::Equal(66),
                CigarOp::Ins(1),
                CigarOp::Equal(116),
                CigarOp::Del(1),
                CigarOp::Equal(25),
                CigarOp::Ins(1),
                CigarOp::Equal(18),
                CigarOp::Ins(1),
                CigarOp::Equal(135),
                CigarOp::Ins(1),
                CigarOp::Equal(88),
                CigarOp::Ins(1),
                CigarOp::Equal(15),
                CigarOp::Ins(1),
                CigarOp::Equal(33),
                CigarOp::Ins(1),
                CigarOp::Equal(139),
                CigarOp::Del(1),
                CigarOp::Equal(11),
                CigarOp::Ins(1),
                CigarOp::Equal(81),
                CigarOp::Ins(1),
                CigarOp::Equal(100),
                CigarOp::Ins(1),
                CigarOp::Equal(131),
                CigarOp::Ins(1),
                CigarOp::Equal(37),
                CigarOp::Ins(1),
                CigarOp::Equal(20),
                CigarOp::Diff(1),
                CigarOp::Equal(94),
                CigarOp::Ins(1),
                CigarOp::Equal(9),
                CigarOp::Del(1),
                CigarOp::Equal(152),
                CigarOp::Ins(1),
                CigarOp::Equal(8),
                CigarOp::Del(1),
                CigarOp::Equal(33),
                CigarOp::Del(1),
                CigarOp::Equal(115),
                CigarOp::Ins(1),
                CigarOp::Equal(198),
                CigarOp::Ins(1),
                CigarOp::Equal(8),
                CigarOp::Ins(1),
                CigarOp::Equal(37),
                CigarOp::Ins(1),
                CigarOp::Equal(14),
                CigarOp::Ins(1),
                CigarOp::Equal(60),
                CigarOp::Diff(1),
                CigarOp::Equal(88),
                CigarOp::Del(1),
                CigarOp::Equal(203),
                CigarOp::Del(1),
                CigarOp::Equal(16),
                CigarOp::Ins(1),
                CigarOp::Equal(24),
                CigarOp::Ins(1),
                CigarOp::Equal(153),
                CigarOp::Del(1),
                CigarOp::Equal(97),
                CigarOp::Ins(1),
                CigarOp::Equal(108),
                CigarOp::Ins(1),
                CigarOp::Equal(9),
                CigarOp::Ins(2),
                CigarOp::Equal(80),
                CigarOp::Del(1),
                CigarOp::Equal(23),
                CigarOp::Ins(1),
                CigarOp::Equal(12),
                CigarOp::Ins(1),
                CigarOp::Equal(171),
                CigarOp::Diff(1),
                CigarOp::Equal(19),
                CigarOp::Del(1),
                CigarOp::Equal(75),
                CigarOp::Ins(1),
                CigarOp::Equal(31),
                CigarOp::Diff(1),
                CigarOp::Equal(3),
                CigarOp::Ins(1),
                CigarOp::Equal(12),
                CigarOp::Ins(1),
                CigarOp::Equal(15),
                CigarOp::Diff(1),
                CigarOp::Equal(47),
                CigarOp::Ins(1),
                CigarOp::Equal(104),
                CigarOp::Ins(2),
                CigarOp::Equal(106),
                CigarOp::Ins(1),
                CigarOp::Equal(4),
                CigarOp::Ins(1),
                CigarOp::Equal(32),
                CigarOp::Del(1),
                CigarOp::Equal(20),
                CigarOp::Ins(1),
                CigarOp::Equal(3),
                CigarOp::Diff(1),
                CigarOp::Equal(61),
                CigarOp::Ins(1),
                CigarOp::Equal(15),
                CigarOp::Diff(1),
                CigarOp::Equal(30),
                CigarOp::Del(1),
                CigarOp::Equal(118),
                CigarOp::Del(1),
                CigarOp::Equal(86),
                CigarOp::Diff(1),
                CigarOp::Equal(28),
                CigarOp::Ins(1),
                CigarOp::Equal(40),
                CigarOp::Diff(1),
                CigarOp::Equal(34),
                CigarOp::Ins(1),
                CigarOp::Equal(25),
                CigarOp::Ins(1),
                CigarOp::Equal(7),
                CigarOp::Ins(1),
                CigarOp::Equal(7),
                CigarOp::Ins(1),
                CigarOp::Equal(36),
                CigarOp::Del(1),
                CigarOp::Equal(100),
                CigarOp::Ins(1),
                CigarOp::Equal(487),
                CigarOp::Ins(1),
                CigarOp::Equal(340),
                CigarOp::Ins(1),
                CigarOp::Equal(198),
                CigarOp::Ins(1),
                CigarOp::Equal(5),
                CigarOp::Ins(1),
                CigarOp::Equal(95),
                CigarOp::Ins(1),
                CigarOp::Equal(203),
                CigarOp::Ins(1),
                CigarOp::Equal(59),
                CigarOp::Ins(1),
                CigarOp::Equal(3),
                CigarOp::Del(1),
                CigarOp::Equal(16),
                CigarOp::Ins(1),
                CigarOp::Equal(30),
                CigarOp::Ins(1),
                CigarOp::Equal(26),
                CigarOp::Ins(1),
                CigarOp::Equal(121),
                CigarOp::Diff(1),
                CigarOp::Equal(114),
                CigarOp::Ins(1),
                CigarOp::Equal(150),
                CigarOp::Ins(1),
                CigarOp::Equal(215),
                CigarOp::Ins(2),
                CigarOp::Equal(209),
                CigarOp::Ins(1),
                CigarOp::Equal(27),
                CigarOp::Ins(1),
                CigarOp::Equal(7),
                CigarOp::Del(1),
                CigarOp::Equal(14),
                CigarOp::Ins(1),
                CigarOp::Equal(410),
                CigarOp::Ins(1),
                CigarOp::Equal(87),
                CigarOp::Ins(1),
                CigarOp::Equal(404),
                CigarOp::Ins(1),
                CigarOp::Equal(19),
                CigarOp::Ins(1),
                CigarOp::Equal(3),
                CigarOp::Ins(1),
                CigarOp::Equal(99),
                CigarOp::Del(1),
                CigarOp::Equal(525),
                CigarOp::Del(27),
                CigarOp::Equal(178),
                CigarOp::Ins(1),
                CigarOp::Equal(177),
                CigarOp::Ins(1),
                CigarOp::Equal(326),
                CigarOp::Diff(1),
                CigarOp::Equal(32),
                CigarOp::Ins(1),
                CigarOp::Equal(13),
                CigarOp::Ins(1),
                CigarOp::Equal(54),
                CigarOp::Del(1),
                CigarOp::Equal(33),
                CigarOp::Ins(1),
                CigarOp::Equal(87),
                CigarOp::Diff(1),
                CigarOp::Equal(39),
                CigarOp::Del(1),
                CigarOp::Equal(90),
                CigarOp::SoftClip(3285),
            ],
        };
        assert_eq!(get_cigar_string(&cigar.ops), CIGAR_STRING);
        assert_eq!(cigar.ops, parse_cigar_string(CIGAR_STRING));
    }

    #[test]
    fn test_mismatch_count_full_cigar() {
        let cigar = Cigar {
            ref_pos: 372,
            ops: parse_cigar_string(CIGAR_STRING),
        };
        let mismatches = extract_snps(
            &cigar,
            &GenomicRegion {
                contig: "chrA".to_string(),
                start: 10002,
                end: 10061,
            },
        );
        assert_eq!(mismatches.len(), 17);
        let ground_truth = vec![
            732, 743, 744, 850, 1567, 3340, 4072, 5061, 5188, 5219, 5537, 5614, 5851, 5920, 7715,
            10709, 10930,
        ];
        assert_eq!(mismatches, ground_truth);
    }

    #[test]
    fn test_mismatch_count_offset_full_cigar() {
        let cigar = Cigar {
            ref_pos: 372,
            ops: parse_cigar_string(CIGAR_STRING),
        };
        let mismatches = extract_snps_offset(
            &cigar,
            &GenomicRegion {
                contig: "chrA".to_string(),
                start: 10002,
                end: 10061,
            },
        );
        assert_eq!(mismatches.len(), 17);
        let ground_truth = vec![
            -9270, -9259, -9258, -9152, -8435, -6662, -5930, -4941, -4814, -4783, -4465, -4388,
            -4151, -4082, -2287, 648, 869,
        ];
        assert_eq!(mismatches, ground_truth);
    }

    #[test]
    fn test_no_mismatch() {
        let cigar = Cigar {
            ref_pos: 0,
            ops: parse_cigar_string("5M"),
        };
        let mismatches = extract_snps(
            &cigar,
            &GenomicRegion {
                contig: "chrA".to_string(),
                start: 10000,
                end: 1000,
            },
        );
        assert_eq!(mismatches.len(), 0);
    }

    #[test]
    fn test_consecutive_mismatches() {
        let cigar = Cigar {
            ref_pos: 0,
            ops: parse_cigar_string("3M4X2M"),
        };
        let mismatches = extract_snps(
            &cigar,
            &GenomicRegion {
                contig: "chrA".to_string(),
                start: 10000,
                end: 1000,
            },
        );
        assert_eq!(mismatches.len(), 4);
        assert_eq!(mismatches, [3, 4, 5, 6]);
    }

    #[test]
    fn test_intermittent_mismatches() {
        let cigar = Cigar {
            ref_pos: 0,
            ops: parse_cigar_string("1M1X1M1X1M1X1M1X1M"),
        };

        let mismatches = extract_snps(
            &cigar,
            &GenomicRegion {
                contig: "chrA".to_string(),
                start: 10000,
                end: 1000,
            },
        );
        assert_eq!(mismatches.len(), 4);
        assert_eq!(mismatches, [1, 3, 5, 7]);

        let mismatches = extract_snps(
            &cigar,
            &GenomicRegion {
                contig: "chrA".to_string(),
                start: 2,
                end: 6,
            },
        );
        assert_eq!(mismatches.len(), 2);
        assert_eq!(mismatches, [1, 7]);
    }
}
