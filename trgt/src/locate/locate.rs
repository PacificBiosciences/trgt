use crate::{reads::HiFiRead, workflows::Params};
use bio::alignment::{pairwise::*, AlignmentOperation};
use itertools::Itertools;
use std::str;

#[derive(Debug, Clone, Copy)]
pub struct TrgtScoring {
    pub match_scr: i32,
    pub mism_scr: i32,
    pub gapo_scr: i32,
    pub gape_scr: i32,
    pub kmer_len: usize,
    pub bandwidth: usize,
}

type Span = (usize, usize);

fn find_spans<F>(
    aligner: &mut banded::Aligner<F>,
    piece: &str,
    seqs: &[&str],
    params: &Params,
) -> Vec<Option<(usize, usize)>>
where
    F: Fn(u8, u8) -> i32,
{
    seqs.iter()
        .map(|s| {
            if let Some(start) = s.find(piece) {
                Some((start, start + piece.len()))
            } else {
                let align = aligner.semiglobal(piece.as_bytes(), s.as_bytes());
                let flank_aln_len = align
                    .operations
                    .iter()
                    .filter(|x| **x == AlignmentOperation::Match)
                    .count();
                if flank_aln_len as f32
                    >= (params.search_flank_len as f32) * params.min_flank_id_frac
                {
                    Some((align.ystart, align.yend))
                } else {
                    None
                }
            }
        })
        .collect()
}

pub fn find_tr_spans(lf: &str, rf: &str, reads: &[HiFiRead], params: &Params) -> Vec<Option<Span>> {
    let lf_piece = &lf[lf.len() - params.search_flank_len..];
    let rf_piece = &rf[..params.search_flank_len];

    let scoring = Scoring {
        match_fn: |a: u8, b: u8| {
            if a == b {
                params.aln_scoring.match_scr
            } else {
                -params.aln_scoring.mism_scr
            }
        },
        match_scores: Some((params.aln_scoring.match_scr, -params.aln_scoring.mism_scr)),
        gap_open: -params.aln_scoring.gapo_scr,
        gap_extend: -params.aln_scoring.gape_scr,
        xclip_prefix: MIN_SCORE,
        xclip_suffix: MIN_SCORE,
        yclip_prefix: 0,
        yclip_suffix: 0,
    };

    let mut aligner = banded::Aligner::with_capacity_and_scoring(
        params.search_flank_len + 10, // global length
        20000,                        // local length: maximum HiFi read length
        scoring,
        params.aln_scoring.kmer_len,
        params.aln_scoring.bandwidth,
    );

    let seqs = reads
        .iter()
        .map(|r| std::str::from_utf8(&r.bases).unwrap())
        .collect_vec();

    let lf_spans = find_spans(&mut aligner, lf_piece, &seqs, params);
    let rf_spans = find_spans(&mut aligner, rf_piece, &seqs, params);

    lf_spans
        .iter()
        .zip(rf_spans.iter())
        .map(|(lf_span, rf_span)| match (lf_span, rf_span) {
            (None, None) => None,      // No left or right span
            (Some(_lf), None) => None, // Left flanking
            (None, Some(_rf)) => None, // Right flanking
            (Some(lf), Some(rf)) => {
                if lf.1 <= rf.0 {
                    Some((lf.1, rf.0))
                } else {
                    None // Discordant flanks
                }
            }
        })
        .collect()
}

/*

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug)]
    struct NoMethRead {
        seq: Seq,
    }

    impl NoMethRead {
        fn new(seq: &str) -> NoMethRead {
            NoMethRead {
                seq: seq.as_bytes().to_vec(),
            }
        }

        fn rc(&self) -> NoMethRead {
            NoMethRead {
                seq: rev_comp(&self.seq),
            }
        }
    }

    impl Read for NoMethRead {
        fn seq(&self) -> &Seq {
            &self.seq
        }

        fn meth_poses(&self) -> Option<&Vec<usize>> {
            None
        }

        fn meth_probs(&self) -> Option<&Vec<u8>> {
            None
        }
    }

    fn get_search_params_for_short_seqs() -> SearchParams {
        SearchParams {
            flank_len: 4,
            kmer_len: 4,
            step_len: 1,
            max_delta: 0,
            min_kmer_count: 2,
        }
    }

    #[test]
    fn extract_spanning_read() {
        let params = get_search_params_for_short_seqs();
        let lf = "AGGCAGGGCG".as_bytes().to_vec();
        let rf = "TCATCGGCTT".as_bytes().to_vec();
        let extractor = RepeatReadExtractor::new(params, lf, rf);

        let query = NoMethRead::new("GGGCGAAAAAAAAAATCATCGG");
        let expected = Some(BasicRead::new("AAAAAAAAAA".as_bytes().to_vec(), None, None));
        assert_eq!(extractor.extract(&query), expected);
        //assert_eq!(extractor.extract(&query.rc()), expected);

        //  //                 -----x          -------
        //  let query = "AGGGCTAAAAAAAAAATCATCGG";
        //  assert_eq!(locator.locate_with_kmers(&query), Some((5, 16, true)));
        //  assert_eq!(locator.locate(&query), Some((6, 16, true)));
        //  assert_eq!(locator.locate_with_kmers(&rc(query)), Some((7, 18, false)));
        //  assert_eq!(locator.locate(&rc(query)), Some((7, 17, false)));

        //  //                 -----x          x------
        //  let query = "AGGGCTAAAAAAAAAAGCATCGG";
        //  assert_eq!(locator.locate_with_kmers(&query), Some((5, 17, true)));
        //  assert_eq!(locator.locate(&query), Some((6, 16, true)));
        //  assert_eq!(locator.locate_with_kmers(&rc(query)), Some((6, 18, false)));
        //  assert_eq!(locator.locate(&rc(query)), Some((7, 17, false)));

        // //                             -----x          x-----
        // let query = "TGGGCTTTTTTTAGGGCTAAAAAAAAAAGCATCG";
        // assert_eq!(locator.locate_with_kmers(&query), Some((17, 29, true)));
        // assert_eq!(locator.locate(&query), Some((18, 28, true)));
        //assert_eq!(locator.locate_with_kmers(&rc(query)), Some((5, 17, false)));
        //assert_eq!(locator.locate(&rc(query)), Some((6, 16, false)));
    }
}


 */
