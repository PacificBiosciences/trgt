use crate::commands::genotype::THREAD_WFA_FLANK;
use crate::trgt::{reads::HiFiRead, workflows::Params};
use rayon::prelude::*;

type Span = (usize, usize);

fn find_spans(piece: &[u8], seqs: &[&[u8]], threshold_frac: f64) -> Vec<Option<(usize, usize)>> {
    seqs.par_iter()
        .map(|s| {
            s.windows(piece.len())
                .position(|window| window == piece)
                .map(|start| (start, start + piece.len()))
                .or_else(|| {
                    THREAD_WFA_FLANK.with(|aligner_cell| {
                        let mut aligner = aligner_cell.borrow_mut();
                        let _status =
                            aligner.align_ends_free(piece, 0, 0, s, s.len() as i32, s.len() as i32);
                        let flank_aln_len = aligner.count_matches() as usize;
                        if flank_aln_len as f64 >= threshold_frac {
                            let ((_pattern_start, _pattern_end), (text_start, text_end)) =
                                aligner.find_alignment_span();
                            Some((text_start, text_end))
                        } else {
                            None
                        }
                    })
                })
        })
        .collect()
}

pub fn find_tr_spans(
    lf: &[u8],
    rf: &[u8],
    reads: &[HiFiRead],
    params: &Params,
) -> Vec<Option<Span>> {
    let lf_piece = &lf[lf.len() - params.search_flank_len..];
    let rf_piece = &rf[..params.search_flank_len];

    let seqs = reads
        .iter()
        .map(|r| r.bases.as_slice())
        .collect::<Vec<&[u8]>>();

    let threshold_frac = (params.search_flank_len as f64) * params.min_flank_id_frac;
    let (lf_spans, rf_spans) = rayon::join(
        || find_spans(lf_piece, &seqs, threshold_frac),
        || find_spans(rf_piece, &seqs, threshold_frac),
    );

    lf_spans
        .into_iter()
        .zip(rf_spans.iter())
        .map(|(lf_span, rf_span)| match (lf_span, rf_span) {
            (None, None) => None,      // No left or right span
            (Some(_lf), None) => None, // Left flanking only
            (None, Some(_rf)) => None, // Right flanking only
            (Some(lf), Some(rf)) => {
                if lf.1 <= rf.0 {
                    Some((lf.1, rf.0))
                } else {
                    None // Discordant flanks (overlap or wrong order)
                }
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_find_vs_windows_comparison() {
        let test_cases = vec![
            // (sequence, piece, expected_span)
            ("ABCDEFG", "CDE", Some((2, 5))),     // Simple match
            ("ABCDEFG", "XYZ", None),             // No match
            ("ABCABCABC", "ABC", Some((0, 3))),   // Multiple matches (finds first)
            ("ABCDEFG", "ABC", Some((0, 3))),     // Match at beginning
            ("ABCDEFG", "EFG", Some((4, 7))),     // Match at end
            ("ABCDEFG", "ABCDEFG", Some((0, 7))), // Piece is the whole sequence
            ("ABC", "ABCDEFG", None),             // Piece longer than sequence
            ("", "ABC", None),                    // Empty sequence
            ("ABCDEFG", "A", Some((0, 1))),       // Single char piece, start
            ("ABCDEFG", "G", Some((6, 7))),       // Single char piece, end
            ("ABCDEFG", "D", Some((3, 4))),       // Single char piece, middle
            ("AAAAA", "AA", Some((0, 2))),        // Overlapping piece in sequence
            ("ACGTNACGT", "N", Some((4, 5))),     // IUPAC codes (valid UTF8)
        ];

        for (seq_str, piece_str, expected_span) in test_cases {
            let seq_bytes = seq_str.as_bytes();
            let piece_bytes = piece_str.as_bytes();

            let find_result = std::str::from_utf8(seq_bytes)
                .unwrap()
                .find(std::str::from_utf8(piece_bytes).unwrap())
                .map(|start| (start, start + piece_bytes.len()));

            let windows_result = if piece_bytes.is_empty() {
                None
            } else {
                seq_bytes
                    .windows(piece_bytes.len())
                    .position(|window| window == piece_bytes)
                    .map(|start| (start, start + piece_bytes.len()))
            };

            assert_eq!(
                find_result, expected_span,
                "Mismatch for find(): seq='{}', piece='{}'",
                seq_str, piece_str
            );

            if !piece_bytes.is_empty() {
                assert_eq!(
                    windows_result, expected_span,
                    "Mismatch for windows(): seq='{}', piece='{}'",
                    seq_str, piece_str
                );
                assert_eq!(
                    find_result, windows_result,
                    "find() and windows() differ: seq='{}', piece='{}'",
                    seq_str, piece_str
                );
            }
        }

        assert_eq!("ABC".find(""), Some(0), "str::find with empty piece");
    }

    // #[derive(Debug)]
    // struct NoMethRead {
    //     seq: Seq,
    // }

    // impl NoMethRead {
    //     fn new(seq: &str) -> NoMethRead {
    //         NoMethRead {
    //             seq: seq.as_bytes().to_vec(),
    //         }
    //     }

    //     fn rc(&self) -> NoMethRead {
    //         NoMethRead {
    //             seq: rev_comp(&self.seq),
    //         }
    //     }
    // }

    // impl Read for NoMethRead {
    //     fn seq(&self) -> &Seq {
    //         &self.seq
    //     }

    //     fn meth_poses(&self) -> Option<&Vec<usize>> {
    //         None
    //     }

    //     fn meth_probs(&self) -> Option<&Vec<u8>> {
    //         None
    //     }
    // }

    // fn get_search_params_for_short_seqs() -> SearchParams {
    //     SearchParams {
    //         flank_len: 4,
    //         kmer_len: 4,
    //         step_len: 1,
    //         max_delta: 0,
    //         min_kmer_count: 2,
    //     }
    // }

    // #[test]
    // fn extract_spanning_read() {
    //     let params = get_search_params_for_short_seqs();
    //     let lf = "AGGCAGGGCG".as_bytes().to_vec();
    //     let rf = "TCATCGGCTT".as_bytes().to_vec();
    //     let extractor = RepeatReadExtractor::new(params, lf, rf);

    //     let query = NoMethRead::new("GGGCGAAAAAAAAAATCATCGG");
    //     let expected = Some(BasicRead::new("AAAAAAAAAA".as_bytes().to_vec(), None, None));
    //     assert_eq!(extractor.extract(&query), expected);
    //     //assert_eq!(extractor.extract(&query.rc()), expected);

    //     //  //                 -----x          -------
    //     //  let query = "AGGGCTAAAAAAAAAATCATCGG";
    //     //  assert_eq!(locator.locate_with_kmers(&query), Some((5, 16, true)));
    //     //  assert_eq!(locator.locate(&query), Some((6, 16, true)));
    //     //  assert_eq!(locator.locate_with_kmers(&rc(query)), Some((7, 18, false)));
    //     //  assert_eq!(locator.locate(&rc(query)), Some((7, 17, false)));

    //     //  //                 -----x          x------
    //     //  let query = "AGGGCTAAAAAAAAAAGCATCGG";
    //     //  assert_eq!(locator.locate_with_kmers(&query), Some((5, 17, true)));
    //     //  assert_eq!(locator.locate(&query), Some((6, 16, true)));
    //     //  assert_eq!(locator.locate_with_kmers(&rc(query)), Some((6, 18, false)));
    //     //  assert_eq!(locator.locate(&rc(query)), Some((7, 17, false)));

    //     // //                             -----x          x-----
    //     // let query = "TGGGCTTTTTTTAGGGCTAAAAAAAAAAGCATCG";
    //     // assert_eq!(locator.locate_with_kmers(&query), Some((17, 29, true)));
    //     // assert_eq!(locator.locate(&query), Some((18, 28, true)));
    //     //assert_eq!(locator.locate_with_kmers(&rc(query)), Some((5, 17, false)));
    //     //assert_eq!(locator.locate(&rc(query)), Some((6, 16, false)));
    // }
}
