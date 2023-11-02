use bio::alignment::{pairwise::*, AlignmentOperation};
use itertools::Itertools;

pub fn repair_consensus(reference: &str, seqs: &[&str], counts: &[usize]) -> String {
    let coverage = counts.iter().sum::<usize>();
    let reference_entry = seqs.iter().zip(counts).find(|rec| *rec.0 == reference);
    let reference_count = match reference_entry {
        Some(rec) => *rec.1,
        None => 0,
    };

    if 2 * reference_count >= coverage {
        return reference.to_string();
    }

    let mut aligner = Aligner::new(-5, -1, |a, b| if a == b { 1i32 } else { -1i32 });
    let aligns = seqs
        .iter()
        .map(|seq| aligner.global(seq.as_bytes(), reference.as_bytes()))
        .collect_vec();

    //                                        A  T  C  G  -
    let mut ref_counts = vec![[0, 0, 0, 0, 0]; reference.len()];
    let mut ref_inserts: Vec<Vec<String>> = vec![Vec::new(); reference.len() + 1];
    for (seq_index, align) in aligns.iter().enumerate() {
        let seq = seqs[seq_index];
        let mut x_pos = 0;
        let mut y_pos = 0;
        for (op, group) in &align.operations.iter().group_by(|a| *a) {
            let op_len = group.count();
            match op {
                AlignmentOperation::Match => {
                    let seq_piece = &seq[x_pos..x_pos + op_len];
                    summarize_matches(&mut ref_counts, y_pos, seq_piece);
                    x_pos += op_len;
                    y_pos += op_len;
                }
                AlignmentOperation::Subst => {
                    let seq_piece = &seq[x_pos..x_pos + op_len];
                    summarize_matches(&mut ref_counts, y_pos, seq_piece);
                    x_pos += op_len;
                    y_pos += op_len;
                }
                AlignmentOperation::Del => {
                    summarize_dels(&mut ref_counts, y_pos, op_len);
                    y_pos += op_len;
                }
                AlignmentOperation::Ins => {
                    let seq_piece = &seq[x_pos..x_pos + op_len];
                    ref_inserts[y_pos].push(seq_piece.to_string());
                    x_pos += op_len;
                }
                _ => panic!("Unexpected operation"),
            };
        }
    }

    let consensus_indexes = ref_counts
        .into_iter()
        .map(|rec| {
            rec.iter()
                .enumerate()
                .max_by_key(|(_, val)| *val)
                .unwrap()
                .0
        })
        .collect_vec();

    let mut consensus: String = String::new();
    for (ref_pos, base_index) in consensus_indexes.iter().enumerate() {
        if ref_inserts[ref_pos].len() > seqs.len() / 2 {
            consensus += get_ins_consensus(&mut ref_inserts[ref_pos], seqs.len());
        }
        if *base_index != 4 {
            consensus.push(match *base_index {
                0 => 'A',
                1 => 'T',
                2 => 'C',
                3 => 'G',
                _ => panic!("Unknown base encoding"),
            });
        }
    }

    consensus
}

fn summarize_matches(ref_counts: &mut [[i32; 5]], ref_pos: usize, seq: &str) {
    for (base_offset, base) in seq.as_bytes().iter().enumerate() {
        let base_index = match base {
            b'A' => 0,
            b'T' => 1,
            b'C' => 2,
            b'G' => 3,
            _ => panic!("Encountered unexpected base"),
        };
        ref_counts[ref_pos + base_offset][base_index] += 1;
    }
}

fn summarize_dels(ref_counts: &mut [[i32; 5]], ref_pos: usize, del_len: usize) {
    let dash_index = 4;
    for base_offset in 0..del_len {
        ref_counts[ref_pos + base_offset][dash_index] += 1;
    }
}

fn get_ins_consensus(ins_by_read: &mut Vec<String>, num_reads: usize) -> &str {
    ins_by_read.sort();
    let reads_without_ins = num_reads - ins_by_read.len();
    let (top_ins, ins_count) = ins_by_read
        .iter()
        .group_by(|ins| *ins)
        .into_iter()
        .map(|(ins, group)| (ins, group.count()))
        .sorted_by(|a, b| Ord::cmp(&b.1, &a.1))
        .next()
        .unwrap();

    if ins_count > reads_without_ins {
        top_ins
    } else {
        ""
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fix_mismatches() {
        let consensus = "CCCCACCCTCCC";
        let sequence1 = "CCCCACCCGCCC";
        let sequence2 = "CCCCCCCCGCCC";
        let sequence3 = "CCCCCCCCCCCC";

        let seqs = vec![sequence1, sequence2, sequence3];
        //let fixed = get_consensus(consensus, &seqs);

        //let expected = "CCCCCCCCGCCC";
        //assert_eq!(fixed, expected);
    }

    #[test]
    fn fix_dels() {
        let consensus = "CCCCCCCCGCCC";
        let sequence1 = "CCCCCCGCCC";
        let sequence2 = "CCCCCCGCCC";
        let sequence3 = "CCCCCCCCGCCC";

        let seqs = vec![sequence1, sequence2, sequence3];
        //let fixed = get_consensus(consensus, &seqs);

        //let expected = "CCCCCCGCCC";
        //assert_eq!(fixed, expected);
    }

    #[test]
    fn fix_ins() {
        let consensus = "CCCCCCCCGCCC";
        let sequence1 = "CCCCCAAACCCGCCC";
        let sequence2 = "CCCCCAAACCCGCCC";
        let sequence3 = "CCCCCACCCGCAACC";

        let seqs = vec![sequence1, sequence2, sequence3];
        //let fixed = crate::locate::get_consensus(consensus, &seqs);

        //let expected = "CCCCCAAACCCGCCC";
        //assert_eq!(fixed, expected);
    }
}
