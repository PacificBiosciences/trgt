use crate::locus::Allele;
use crate::locus::Locus;
use crate::read::Read;
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;

#[derive(Debug, Clone)]
pub struct AlignInfo {
    pub seq: String,
    pub meth: Option<Vec<u8>>,
    pub align: Alignment,
    pub allele_index: usize,
}

pub fn align_reads(genotype: &Vec<Allele>, reads: Vec<Read>) -> Vec<AlignInfo> {
    let likely_max_len = get_longest_allele(genotype) + 50;
    let mut aligner = get_aligner(likely_max_len);
    let mut align_infos = Vec::new();

    for read in reads {
        let i = read.allele as usize;
        align_infos.push(AlignInfo {
            align: aligner.global(read.seq.as_bytes(), genotype[i].seq.as_bytes()),
            seq: read.seq,
            meth: read.meth,
            allele_index: i,
        });
    }

    align_infos
}

pub fn align_to_flanks(locus: &Locus, reads: Vec<Read>) -> Vec<AlignInfo> {
    let max_read_len = reads.iter().map(|r| r.seq.len()).max().unwrap();
    let flank_len = locus.left_flank.len();

    let mut aligns = Vec::new();
    let mut aligner = get_aligner(max_read_len);
    for read in reads {
        let mut ref_seq = locus.left_flank[locus.left_flank.len() - flank_len..].to_string();
        ref_seq += &read.seq[read.left_flank..read.seq.len() - read.right_flank];
        ref_seq += &locus.right_flank[..flank_len];

        let align = aligner.global(read.seq.as_bytes(), ref_seq.as_bytes());

        aligns.push(AlignInfo {
            seq: read.seq,
            meth: read.meth,
            align,
            allele_index: 0,
        });
    }
    aligns
}

type ScoreFunc = fn(u8, u8) -> i32;

fn score(a: u8, b: u8) -> i32 {
    if a == b {
        1i32
    } else {
        -1i32
    }
}

fn get_aligner<'a>(likely_max_len: usize) -> Aligner<&'a ScoreFunc> {
    Aligner::with_capacity(
        likely_max_len,
        likely_max_len,
        -5,
        -1,
        &(score as ScoreFunc),
    )
}

fn get_longest_allele(genotype: &Vec<Allele>) -> usize {
    let mut max_len = 0;

    for allele in genotype {
        if allele.seq.len() > max_len {
            max_len = allele.seq.len();
        }
    }

    max_len
}
