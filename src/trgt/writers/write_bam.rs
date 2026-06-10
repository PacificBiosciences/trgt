//! Defines the `BamWriter` struct and associated functions for creating and writing spanning reads to a BAM file.
//!
use crate::cli;
use crate::trgt::{
    locus::Locus, reads::extract_mismatch_offsets, reads::AlleleAssign, workflows::LocusResult,
};
use crate::utils::Result;
use rust_htslib::bam::{
    self,
    header::HeaderRecord,
    record::{Aux, CigarString},
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
        n_compression_threads: usize,
    ) -> Result<BamWriter> {
        let header = Self::create_header(template_header);
        let mut writer = bam::Writer::from_path(output_bam_path, &header, bam::Format::Bam)
            .map_err(|e| e.to_string())?;
        writer.set_threads(n_compression_threads).unwrap();
        log::debug!("Creating BAM writer");
        Ok(BamWriter { writer, flank_len })
    }

    /// Creates a BAM header based on a template, including additional program information.
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

    #[inline]
    fn assign_to_i32(assign: AlleleAssign) -> i32 {
        match assign {
            AlleleAssign::A0 => 0,
            AlleleAssign::A1 => 1,
            AlleleAssign::Both => 2,
            AlleleAssign::None => -1,
        }
    }

    #[inline]
    fn contig_tid(&self, contig: &str) -> i32 {
        self.writer
            .header()
            .tid(contig.as_bytes())
            .expect("unknown contig") as i32
    }

    /// Writes the spanning reads to the BAM file from the genotyping results for a specific locus.
    ///
    /// # Arguments
    ///
    /// * `locus` - `Locus` struct containing locus information.
    /// * `results` - `LocusResult` struct containing genotyping results.
    pub fn write(&mut self, locus: &Locus, results: &LocusResult) {
        if results.spanning_reads.is_empty() {
            return;
        }

        let tid = self.contig_tid(&locus.region.contig);
        for sr in &results.spanning_reads {
            let read = &sr.read;
            let span = &sr.span;

            if span.0 < self.flank_len || read.bases.len() < span.1 + self.flank_len {
                log::error!("Read {} has unexpectedly short flanks", read.id);
                continue;
            }

            let left_clip_len = span.0 - self.flank_len;
            let right_clip_len = read.bases.len() - span.1 - self.flank_len;

            let read = match read.clip_bases(left_clip_len, right_clip_len) {
                Some(r) => r,
                None => {
                    log::error!("Read {} has unexpectedly short flanks", read.id);
                    continue;
                }
            };

            let mut rec = bam::Record::new();
            rec.set_tid(tid);

            if read.is_reverse {
                rec.set_reverse();
            }

            if let Some(cigar) = read.cigar {
                rec.set_pos(cigar.ref_pos);
                rec.set(
                    read.id.as_bytes(),
                    Some(&CigarString(cigar.ops)),
                    &read.bases,
                    &read.quals,
                );
                rec.set_mapq(read.mapq);
                rec.unset_unmapped();
            } else {
                rec.set_pos(locus.region.start as i64);
                rec.set(read.id.as_bytes(), None, &read.bases, &read.quals);
                rec.set_unmapped();
            }

            rec.push_aux(b"TR", Aux::String(&locus.id)).unwrap();
            rec.push_aux(b"rq", Aux::Float(read.read_qual.unwrap_or(-1.0f32)))
                .unwrap();

            if let Some(meth) = &read.meth {
                rec.push_aux(b"MC", Aux::ArrayU8(meth.into())).unwrap();
            }

            let mismatch_offsets = read.mismatch_positions.as_ref().map(|positions| {
                extract_mismatch_offsets(positions, locus.region.start, locus.region.end)
            });

            if let Some(mismatches) = mismatch_offsets {
                rec.push_aux(b"MO", Aux::ArrayI32((&mismatches).into()))
                    .unwrap();
            }

            if let Some(hp) = read.hp_tag {
                rec.push_aux(b"HP", Aux::U8(hp)).unwrap();
            }

            if let Some(ps) = read.ps_tag {
                rec.push_aux(b"PS", Aux::U32(ps)).unwrap();
            }

            let start_offset = (read.ref_start - locus.region.start as i64) as i32;
            let end_offset = (read.ref_end - locus.region.end as i64) as i32;

            rec.push_aux(b"SO", Aux::I32(start_offset)).unwrap();
            rec.push_aux(b"EO", Aux::I32(end_offset)).unwrap();

            rec.push_aux(b"AL", Aux::I32(Self::assign_to_i32(sr.assign)))
                .unwrap();

            rec.push_aux(
                b"FL",
                Aux::ArrayU32((&[self.flank_len as u32, self.flank_len as u32]).into()),
            )
            .unwrap();

            self.writer.write(&rec).unwrap();
        }
    }
}
