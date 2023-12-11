use crate::cli;
use crate::locus::Locus;
use crate::reads::clip_bases;
use crate::workflows::LocusResult;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{AuxArray, CigarString};
use rust_htslib::{bam, bam::record::Aux};
use std::env;

pub struct BamWriter {
    writer: bam::Writer,
    output_flank_len: usize,
}

impl BamWriter {
    fn update_header(header: bam::Header) -> bam::Header {
        let mut header = header;
        let args: Vec<String> = env::args().collect();
        let command_line = args.join(" ");

        let mut record = HeaderRecord::new(b"PG");
        record.push_tag(b"ID", env!("CARGO_PKG_NAME"));
        record.push_tag(b"PN", env!("CARGO_PKG_NAME"));
        record.push_tag(b"CL", command_line);
        record.push_tag(b"VN", (*cli::FULL_VERSION).to_string());
        header.push_record(&record);

        header
    }

    pub fn new(
        output_bam_path: &str,
        bam_header: bam::Header,
        output_flank_len: usize,
    ) -> Result<BamWriter, String> {
        let bam_header = Self::update_header(bam_header);
        let writer = bam::Writer::from_path(output_bam_path, &bam_header, bam::Format::Bam)
            .map_err(|e| e.to_string())?;
        Ok(BamWriter {
            writer,
            output_flank_len,
        })
    }

    pub fn write(&mut self, locus: &Locus, results: &LocusResult) {
        let num_reads = results.reads.len();
        for index in 0..num_reads {
            let read = &results.reads[index];
            let classification = results.classification[index];
            let span = &results.tr_spans[index];

            if span.0 < self.output_flank_len || read.bases.len() < span.1 + self.output_flank_len {
                log::error!("Read {} has unexpectedly short flanks", read.id);
                continue;
            }

            let left_clip_len = span.0 - self.output_flank_len;
            let right_clip_len = read.bases.len() - span.1 - self.output_flank_len;
            let clipped_read = clip_bases(read, left_clip_len, right_clip_len);
            if clipped_read.is_none() {
                log::error!("Read {} has unexpectedly short flanks", read.id);
                continue;
            }
            let read = clipped_read.unwrap();

            let contig = locus.region.contig.as_bytes();
            let contig_id = self.writer.header().tid(contig).unwrap();
            let quals = "(".repeat(read.bases.len());

            let mut rec = bam::Record::new();
            rec.set_tid(contig_id as i32);

            if let Some(cigar) = read.cigar {
                rec.set_pos(cigar.ref_pos);
                let cigar = CigarString(cigar.ops);
                rec.set(
                    read.id.as_bytes(),
                    Some(&cigar),
                    &read.bases,
                    quals.as_bytes(),
                );
                rec.set_mapq(read.mapq);
            } else {
                rec.set(read.id.as_bytes(), None, &read.bases, quals.as_bytes());
                rec.set_pos(locus.region.start as i64);
                rec.set_unmapped();
            }

            let tr_tag = Aux::String(&locus.id);
            rec.push_aux(b"TR", tr_tag).unwrap();

            if let Some(meth) = &read.meth {
                let mc_tag = Aux::ArrayU8(meth.into());
                rec.push_aux(b"MC", mc_tag).unwrap();
            }

            if let Some(mismatch) = &read.mismatch_offsets {
                let mm_tag = Aux::ArrayI32(mismatch.into());
                rec.push_aux(b"MO", mm_tag).unwrap();
            }

            if let Some(hp) = read.hp_tag {
                let hp_tag = Aux::U8(hp);
                rec.push_aux(b"HP", hp_tag).unwrap();
            }

            rec.push_aux(b"SO", Aux::I32(read.start_offset)).unwrap();
            rec.push_aux(b"EO", Aux::I32(read.end_offset)).unwrap();
            rec.push_aux(b"AL", Aux::I32(classification)).unwrap();

            let dat: &Vec<u32> = &vec![self.output_flank_len as u32, self.output_flank_len as u32];
            let fl_tag: AuxArray<u32> = dat.into();
            rec.push_aux(b"FL", Aux::ArrayU32(fl_tag)).unwrap();

            self.writer.write(&rec).unwrap();
        }
    }
}
