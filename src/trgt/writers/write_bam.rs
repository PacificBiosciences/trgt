//! Defines the `BamWriter` struct and associated functions for creating and writing spanning reads to a BAM file.
//!
use crate::cli;
use crate::trgt::{locus::Locus, workflows::LocusResult};
use crate::utils::Result;
use rust_htslib::bam::{
    self,
    header::HeaderRecord,
    record::{Aux, AuxArray, CigarString},
};
use std::env;

/// Structure for writing spanning reads from genotyping results.
pub struct BamWriter {
    /// The BAM writer used to write reads to the BAM file.
    writer: bam::Writer,
    /// The length of the flanking sequence to include on each side of the read when writing to the BAM file.
    flank_len: usize,
}

impl BamWriter {
    /// Constructs a new `BamWriter` instance.
    ///
    /// # Arguments
    /// * `output_bam_path` - Path of the output BAM file.
    /// * `template_header` - The BAM header to use as a template for the output file.
    /// * `flank_len` - The length of the flanking sequence to include on each side of the read.
    ///
    /// # Returns
    /// Returns a `Result` with either a new `BamWriter` instance or an error message.
    pub fn new(
        output_bam_path: &str,
        template_header: bam::Header,
        flank_len: usize,
    ) -> Result<BamWriter> {
        let header = Self::create_header(template_header);
        let writer = bam::Writer::from_path(output_bam_path, &header, bam::Format::Bam)
            .map_err(|e| e.to_string())?;
        Ok(BamWriter { writer, flank_len })
    }

    /// Creates a BAM header based on a template, including additional program information..
    ///
    /// # Arguments
    ///
    /// * `template_header` - BAM header to be used as template.
    ///
    /// # Returns
    ///
    /// Returns the updated BAM header.
    fn create_header(template_header: bam::Header) -> bam::Header {
        let mut header = template_header;
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

    /// Writes the spanning reads to the BAM file from the genotyping results for a specific locus.
    ///
    /// # Arguments
    ///
    /// * `locus` - `Locus` struct containing locus information.
    /// * `results` - `LocusResult` struct containing genotyping results.
    pub fn write(&mut self, locus: &Locus, results: &LocusResult) {
        let num_reads = results.reads.len();
        for index in 0..num_reads {
            let read = &results.reads[index];
            let classification = results.classification[index];
            let span = &results.tr_spans[index];

            if span.0 < self.flank_len || read.bases.len() < span.1 + self.flank_len {
                log::error!("Read {} has unexpectedly short flanks", read.id);
                continue;
            }

            let left_clip_len = span.0 - self.flank_len;
            let right_clip_len = read.bases.len() - span.1 - self.flank_len;
            let clipped_read = read.clip_bases(left_clip_len, right_clip_len);
            if clipped_read.is_none() {
                log::error!("Read {} has unexpectedly short flanks", read.id);
                continue;
            }
            let read = clipped_read.unwrap();

            let contig = locus.region.contig.as_bytes();
            let contig_id = self.writer.header().tid(contig).unwrap();

            let mut rec = bam::Record::new();
            rec.set_tid(contig_id as i32);
            if read.is_reverse {
                rec.set_reverse();
            }

            if let Some(cigar) = read.cigar {
                rec.set_pos(cigar.ref_pos);
                let cigar = CigarString(cigar.ops);
                rec.set(read.id.as_bytes(), Some(&cigar), &read.bases, &read.quals);
                rec.set_mapq(read.mapq);
            } else {
                rec.set(read.id.as_bytes(), None, &read.bases, &read.quals);
                rec.set_pos(locus.region.start as i64);
                rec.set_unmapped();
            }

            let tr_tag = Aux::String(&locus.id);
            rec.push_aux(b"TR", tr_tag).unwrap();

            rec.push_aux(b"rq", Aux::Float(read.read_qual.unwrap_or(-1.0) as f32))
                .unwrap();

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

            let dat: &Vec<u32> = &vec![self.flank_len as u32, self.flank_len as u32];
            let fl_tag: AuxArray<u32> = dat.into();
            rec.push_aux(b"FL", Aux::ArrayU32(fl_tag)).unwrap();

            self.writer.write(&rec).unwrap();
        }
    }
}
