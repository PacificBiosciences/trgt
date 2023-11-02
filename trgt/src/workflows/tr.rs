use crate::cluster;
use crate::genotype::{self, flank_genotype, Gt};
use crate::label::label_alleles;
use crate::locate::{Locator, TrgtScoring};
use crate::locus::{Genotyper, Locus};
use crate::reads::{clip_to_region, HiFiRead};
use crate::snp;
use crate::workflows::{Allele, Genotype, LocusResult};
use itertools::Itertools;
use rust_htslib::bam;
use std::vec;

pub type Result<T> = std::result::Result<T, String>;

pub struct Params {
    pub search_flank_len: usize,
    pub min_read_qual: f64,
    pub max_depth: usize,
    pub aln_scoring: TrgtScoring,
    pub min_flank_id_frac: f32,
}

pub fn analyze(
    locus: &Locus,
    params: &Params,
    bam: &mut bam::IndexedReader,
) -> Result<LocusResult> {
    if locus.ploidy == genotype::Ploidy::Zero {
        return Ok(LocusResult::empty());
    }
    let mut reads = extract_reads(
        locus,
        params.search_flank_len as u32,
        bam,
        params.min_read_qual,
    )?;
    log::debug!("{}: Collected {} reads", locus.id, reads.len());

    snp::analyze_snps(&mut reads, &locus.region);

    let clip_radius = 500;
    let reads = clip_reads(locus, clip_radius, reads);
    log::debug!("{}: {} reads left after clipping", locus.id, reads.len());

    let (reads, spans) = get_spanning_reads(locus, params, reads);
    if reads.is_empty() {
        return Ok(LocusResult::empty());
    }

    let trs = reads
        .iter()
        .zip(spans.iter())
        .map(|(r, s)| std::str::from_utf8(&r.bases[s.0..s.1]).unwrap())
        .collect_vec();

    let (mut gt, mut allele_seqs, mut classification) = match locus.genotyper {
        Genotyper::Size => genotype::genotype(locus.ploidy, &trs),
        Genotyper::Cluster => {
            let spanning = reads.iter().map(|r| &r.bases[..]).collect_vec();
            cluster::genotype(&spanning, &trs)
        }
    };

    // Attempt flank re-genotyping only if alleles have similar length
    if gt.len() == 2 && gt[0].size.abs_diff(gt[1].size) <= 10 {
        let snp_result = flank_genotype(&reads, &trs);
        if let Some((snp_gt, snp_alleles, snp_assignment)) = snp_result {
            (gt, allele_seqs, classification) = (snp_gt, snp_alleles, snp_assignment);
        }
    }

    let annotations = label_alleles(locus, &allele_seqs);

    let spanning_by_hap = [
        classification.iter().filter(|&x| *x == 0).count(),
        classification.iter().filter(|&x| *x == 1).count(),
    ];
    let meth_by_hap = get_meth(&gt, &reads, &spans);
    let mut genotype = Genotype::new();
    for allele_index in 0..gt.len() {
        genotype.push(Allele {
            seq: allele_seqs[allele_index].clone(),
            annotation: annotations[allele_index].clone(),
            ci: gt[allele_index].ci,
            num_spanning: spanning_by_hap[allele_index],
            meth: meth_by_hap[allele_index],
        });
    }

    // Put reference allele first
    if genotype.len() != 1 && genotype[0].seq != locus.tr && genotype[1].seq == locus.tr {
        genotype.swap(0, 1);
        for c in classification.iter_mut() {
            *c = 1 - *c;
        }
    }

    Ok(LocusResult {
        genotype,
        reads,
        tr_spans: spans,
        classification,
    })
}

fn get_spanning_reads(
    locus: &Locus,
    params: &Params,
    reads: Vec<HiFiRead>,
) -> (Vec<HiFiRead>, Vec<(usize, usize)>) {
    let mut locator = Locator::new(
        &locus.left_flank,
        &locus.right_flank,
        params.search_flank_len,
    );
    let tr_spans = locator.locate(&reads, params);

    let reads_and_spans = reads
        .into_iter()
        .zip(tr_spans)
        .filter(|(_r, s)| s.is_some())
        .map(|(r, s)| (r, s.unwrap()))
        .collect_vec();
    log::debug!(
        "{}: Found {} spanning reads",
        locus.id,
        reads_and_spans.len()
    );

    // Check flanks
    let mut reads_and_spans = reads_and_spans
        .into_iter()
        .filter(|(r, s)| {
            s.0 >= params.search_flank_len && r.bases.len() - s.1 >= params.search_flank_len
        })
        .collect_vec();
    log::debug!(
        "{}: {} spanning reads had sufficiently long flanks",
        locus.id,
        reads_and_spans.len()
    );

    if reads_and_spans.is_empty() {
        return (Vec::new(), Vec::new());
    }

    // Sort and downsample
    reads_and_spans.sort_by(|(_ra, sa), (_rb, sb)| (sa.1 - sa.0).cmp(&(sb.1 - sb.0)));
    uniform_downsample(&mut reads_and_spans, params.max_depth);
    log::debug!(
        "{}: downsampled to {} reads",
        locus.id,
        reads_and_spans.len()
    );

    let (reads, spans): (Vec<_>, Vec<_>) = reads_and_spans.into_iter().unzip();

    (reads, spans)
}

fn uniform_downsample(reads_and_spans: &mut Vec<(HiFiRead, (usize, usize))>, output_length: usize) {
    let num_reads = reads_and_spans.len();
    if num_reads > output_length {
        let mut fast: f64 = 0.0;
        let step = (num_reads as f64) / (output_length as f64);
        for i in 0..output_length {
            let ind = fast.floor() as usize;
            if ind != i {
                reads_and_spans.swap(i, ind);
            }
            fast += step;
        }
        reads_and_spans.truncate(output_length);
    }
}

fn clip_reads(locus: &Locus, radius: usize, reads: Vec<HiFiRead>) -> Vec<HiFiRead> {
    let region = (
        locus.region.start as i64 - radius as i64,
        locus.region.end as i64 + radius as i64,
    );

    reads
        .into_iter()
        .filter_map(|r| clip_to_region(r, region))
        .collect_vec()
}

fn get_meth(gt: &Gt, reads: &[HiFiRead], spans: &[(usize, usize)]) -> Vec<Option<f64>> {
    let mut meths_1 = Vec::new();
    let mut meths_2 = Vec::new();

    for (read, span) in reads.iter().zip(spans.iter()) {
        if read.meth.is_none() {
            continue;
        }
        let level = get_tr_meth(read, span);
        if level.is_none() {
            continue;
        }
        let level = level.unwrap();
        let assignment = assign_read(gt, span.1 - span.0);
        match assignment {
            Assignment::First => meths_1.push(level),
            Assignment::Second => meths_2.push(level),
            Assignment::Both => {
                meths_1.push(level);
                meths_2.push(level);
            }
            Assignment::None => {}
        };
    }

    let meth_1 = if !meths_1.is_empty() {
        let meth = meths_1.iter().sum::<f64>() / meths_1.len() as f64;
        Some(meth)
    } else {
        None
    };
    let meth_2 = if !meths_2.is_empty() {
        let meth = meths_2.iter().sum::<f64>() / meths_2.len() as f64;
        Some(meth)
    } else {
        None
    };

    if gt.len() == 2 {
        vec![meth_1, meth_2]
    } else {
        vec![meth_1]
    }
}

enum Assignment {
    First,
    Second,
    Both,
    None,
}

fn assign_read(gt: &Gt, tr_len: usize) -> Assignment {
    if gt.len() == 1 {
        return Assignment::First;
    }

    let hap1_len = gt[0].size;
    let hap2_len = gt[1].size;

    let spans_1 = gt[0].ci.0 <= tr_len && tr_len <= gt[0].ci.1;
    let spans_2 = gt[1].ci.0 <= tr_len && tr_len <= gt[1].ci.1;

    let dist_1 = tr_len.abs_diff(gt[0].size);
    let dist_2 = tr_len.abs_diff(gt[1].size);

    if dist_1 < dist_2 && spans_1 {
        return Assignment::First;
    }

    if dist_2 < dist_1 && spans_2 {
        return Assignment::Second;
    }

    if hap1_len == hap2_len && spans_1 {
        return Assignment::Both;
    }

    Assignment::None
}

fn extract_reads(
    locus: &Locus,
    flank_len: u32,
    bam: &mut bam::IndexedReader,
    min_read_qual: f64,
) -> Result<Vec<HiFiRead>> {
    let mut reads = Vec::new();
    let extraction_region = (
        locus.region.contig.as_str(),
        locus.region.start.saturating_sub(flank_len),
        locus.region.end + flank_len,
    );
    if let Err(msg) = bam.fetch(extraction_region) {
        log::warn!("{}", msg);
        return Ok(Vec::new());
    }

    let mut num_filtered = 0;
    for rec in bam::Read::records(bam) {
        let rec = rec.map_err(|e| e.to_string())?;
        if rec.is_supplementary() || rec.is_secondary() {
            continue;
        }

        let read = HiFiRead::from_hts_rec(rec, &locus.region);
        if let Some(qual) = read.read_qual {
            if qual >= min_read_qual {
                reads.push(read);
            } else {
                num_filtered += 1;
            }
        } else {
            reads.push(read);
        }
    }

    if num_filtered > 0 {
        let total = num_filtered + reads.len();
        log::warn!("Quality filtered {} out of {} reads", num_filtered, total);
    }
    Ok(reads)
}

fn get_tr_meth(read: &HiFiRead, span: &(usize, usize)) -> Option<f64> {
    if read.meth.is_none() || read.meth.as_ref().unwrap().is_empty() {
        return None;
    }

    let meth = read.meth.as_ref().unwrap();

    let mut total_meth = 0.0;
    let mut cpg_count = 0;
    let mut cpg_index = 0;
    for pos in 0..read.bases.len() - 1 {
        let dinuc = &read.bases[pos..pos + 2];

        if dinuc == b"CG" {
            if span.0 <= pos && pos < span.1 {
                cpg_count += 1;
                total_meth += match meth.get(cpg_index) {
                    Some(value) => *value as f64 / 255.0,
                    None => {
                        log::error!("Read {} has malformed methylation profile", read.id);
                        std::process::exit(1);
                    }
                };
            }

            cpg_index += 1;
        }
    }

    if cpg_count != 0 {
        let mean_meth = total_meth / cpg_count as f64;
        Some(mean_meth)
    } else {
        None
    }
}
