use super::{
    strategy::exact::merge_exact,
    vcf_reader::{VcfReader, VcfReaders},
    vcf_writer::VcfWriter,
};
use crate::{
    cli::MergeArgs,
    utils::{format_number_with_commas, open_genome_reader, Result},
};
use once_cell::sync::Lazy;
use rust_htslib::{
    bcf::{self, record::GenotypeAllele, Record},
    faidx,
};
use semver::Version;
use std::{
    any::Any,
    cmp::Ordering,
    collections::{BinaryHeap, HashSet},
    env,
    path::PathBuf,
};

const MISSING_INTEGER: i32 = i32::MIN;
const VECTOR_END_INTEGER: i32 = i32::MIN + 1;
static MISSING_FLOAT: Lazy<f32> = Lazy::new(|| f32::from_bits(0x7F80_0001));
static VECTOR_END_FLOAT: Lazy<f32> = Lazy::new(|| f32::from_bits(0x7F80_0002));

struct FormatData {
    als: Vec<i32>,
    allrs: Vec<Vec<u8>>,
    sds: Vec<i32>,
    mcs: Vec<Vec<u8>>,
    mss: Vec<Vec<u8>>,
    aps: Vec<f32>,
    ams: Vec<f32>,
}

impl FormatData {
    fn new() -> Self {
        FormatData {
            als: Vec::new(),
            allrs: Vec::new(),
            sds: Vec::new(),
            mcs: Vec::new(),
            mss: Vec::new(),
            aps: Vec::new(),
            ams: Vec::new(),
        }
    }
}

trait PushMissingAndEnd: Any {
    fn missing() -> Self;
    fn vector_end() -> Self;

    fn push_missing_and_end(vec: &mut Vec<Self>)
    where
        Self: Sized,
    {
        vec.push(Self::missing());
        vec.push(Self::vector_end());
    }
}

macro_rules! impl_push_missing_and_end {
    ($type:ty, $missing:expr, $end:expr) => {
        impl PushMissingAndEnd for $type {
            fn missing() -> Self {
                $missing
            }

            fn vector_end() -> Self {
                $end
            }
        }
    };
    ($type:ty, $missing:expr, $end:expr, $custom_push:expr) => {
        impl PushMissingAndEnd for $type {
            fn missing() -> Self {
                $missing
            }

            fn vector_end() -> Self {
                $end
            }

            #[allow(clippy::redundant_closure_call)]
            fn push_missing_and_end(vec: &mut Vec<Self>) {
                ($custom_push)(vec);
            }
        }
    };
}

impl_push_missing_and_end!(i32, MISSING_INTEGER, VECTOR_END_INTEGER);
impl_push_missing_and_end!(f32, *MISSING_FLOAT, *VECTOR_END_FLOAT);
impl_push_missing_and_end!(Vec<u8>, Vec::new(), Vec::new(), |vec: &mut Vec<Vec<u8>>| {
    vec.push(vec![b'.']);
});

enum FieldType {
    String,
    Integer,
}

#[derive(Debug)]
struct VcfRecordWithSource {
    record: bcf::Record,
    reader_index: usize,
}

impl PartialEq for VcfRecordWithSource {
    fn eq(&self, other: &Self) -> bool {
        self.record.pos() == other.record.pos()
    }
}

impl Eq for VcfRecordWithSource {}

impl PartialOrd for VcfRecordWithSource {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for VcfRecordWithSource {
    fn cmp(&self, other: &Self) -> Ordering {
        self.record.pos().cmp(&other.record.pos()).reverse()
    }
}

pub struct VcfProcessor {
    pub readers: Vec<VcfReader>,
    pub writer: VcfWriter,
    pub genome_reader: Option<faidx::Reader>,
    pub contig_order: Vec<String>,
    pub skip_n: usize,
    pub process_n: usize,
    pub needs_padding: bool,
    pub quit_on_error: bool,
}

impl VcfProcessor {
    pub fn new(args: &MergeArgs, vcfs: Vec<PathBuf>) -> Result<Self> {
        let vcf_readers = VcfReaders::new(vcfs)?;
        if vcf_readers.readers.len() == 1 && !args.force_single {
            return Err("Expected two or more files to merge, got only one. Use --force-single to proceed anyway".into());
        }

        let mut contig_order = vcf_readers.get_contig_order()?;

        if let Some(ref user_contigs) = args.contigs {
            let user_contig_set: HashSet<&String> = user_contigs.iter().collect();
            let original_contig_set: HashSet<&String> = contig_order.iter().collect();

            if !user_contig_set.is_subset(&original_contig_set) {
                let missing_contigs: Vec<&&String> =
                    user_contig_set.difference(&original_contig_set).collect();
                return Err(format!(
                    "The following user-specified contigs do not exist in the VCF files: {:?}",
                    missing_contigs
                ));
            }
            contig_order.retain(|contig| user_contig_set.contains(contig));
        }

        let needs_padding = vcf_readers
            .readers
            .iter()
            .any(|reader| reader.version.major < Version::new(1, 0, 0).major);

        let genome_reader = if needs_padding {
            Some(open_genome_reader(args.genome_path.as_ref().ok_or(
                "A reference genome is required for merging pre v1.0 TRGT VCFs, provide as --genome ref.fa"
            )?)?)
        } else {
            None
        };

        let out_header = Self::create_output_header(&vcf_readers, args)?;
        let writer = VcfWriter::new(&out_header, &args.output_type, args.output.as_ref())?;

        if needs_padding {
            log::debug!("At least one VCF file is pre-1.0 and needs base padding!");
        }

        Ok(VcfProcessor {
            readers: vcf_readers.readers,
            writer,
            genome_reader,
            contig_order,
            skip_n: args.skip_n.unwrap_or(0),
            process_n: args.process_n.unwrap_or(usize::MAX),
            needs_padding,
            quit_on_error: args.quit_on_error,
        })
    }

    fn create_output_header(vcf_readers: &VcfReaders, args: &MergeArgs) -> Result<bcf::Header> {
        let mut out_header = bcf::Header::new();
        vcf_readers.merge_headers(&mut out_header, args.force_samples)?;

        // Update header fields to be consistent
        out_header.remove_format(b"ALCI");
        out_header.remove_format(b"AM");
        out_header.push_record(
            b"##FORMAT=<ID=ALLR,Number=.,Type=String,Description=\"Length range per allele\">",
        );
        out_header.push_record(
            b"##FORMAT=<ID=AM,Number=.,Type=Float,Description=\"Mean methylation level per allele\">",
        );

        if !args.no_version {
            Self::add_version_info(&mut out_header);
        }

        Ok(out_header)
    }

    fn add_version_info(out_header: &mut bcf::Header) {
        let version_line = format!(
            "##{}Version={}",
            env!("CARGO_PKG_NAME"),
            *crate::cli::FULL_VERSION
        );
        out_header.push_record(version_line.as_bytes());

        let command_line = format!(
            "##{}Command={}",
            env!("CARGO_PKG_NAME"),
            env::args().collect::<Vec<String>>().join(" ")
        );
        out_header.push_record(command_line.as_bytes());
    }

    pub fn merge_variants(&mut self) -> Result<()> {
        let mut n = 0;
        let mut n_processed = 0;
        let mut n_failed = 0;

        let mut sample_records = vec![None; self.readers.len()];

        for contig in &self.contig_order.clone() {
            let mut heap = self.init_heap(contig)?;

            while let Some(min_element) = heap.pop() {
                let min_pos = min_element.record.pos();
                sample_records[min_element.reader_index] = Some(min_element.record);

                while let Some(peek_next_element) = heap.peek() {
                    if peek_next_element.record.pos() == min_pos {
                        let next_element = heap.pop().unwrap();
                        sample_records[next_element.reader_index] = Some(next_element.record);
                    } else {
                        break;
                    }
                }

                if n >= self.skip_n {
                    log::trace!("Processing: {}:{}", contig, min_pos);
                    if self.needs_padding {
                        self.add_padding_base(&mut sample_records, contig, min_pos);
                    }
                    match self.merge_variant(&sample_records, contig, min_pos) {
                        Ok(_) => {
                            n_processed += 1;
                            if n_processed >= self.process_n {
                                return Ok(());
                            }
                        }
                        Err(e) => {
                            if self.quit_on_error {
                                return Err(e);
                            }
                            n_failed += 1;
                            log::warn!("{} Skipping...", e);
                        }
                    }
                }
                n += 1;

                self.update_heap(&mut heap, &sample_records)?;
                sample_records.fill(None);
            }
        }
        let mut log_message = format!(
            "Successfully merged {} TR sites.",
            format_number_with_commas(n_processed)
        );
        if n_failed > 0 {
            log_message.push_str(&format!(
                " Failed to merge {} TR sites!",
                format_number_with_commas(n_failed)
            ));
        }
        log::info!("{}", log_message);
        Ok(())
    }

    fn init_heap(&mut self, contig: &str) -> Result<BinaryHeap<VcfRecordWithSource>> {
        let mut heap = BinaryHeap::new();
        for (index, reader) in self.readers.iter_mut().enumerate() {
            let rid = match reader.header.name2rid(contig.as_bytes()) {
                Ok(id) => id,
                Err(_) => continue,
            };

            if reader.reader.fetch(rid, 0, None).is_err() {
                continue;
            }

            if reader.advance() {
                heap.push(VcfRecordWithSource {
                    record: reader.current_record.clone(),
                    reader_index: index,
                });
            }
        }
        Ok(heap)
    }

    fn update_heap(
        &mut self,
        heap: &mut BinaryHeap<VcfRecordWithSource>,
        sample_records: &[Option<Record>],
    ) -> Result<()> {
        for (index, record) in sample_records.iter().enumerate() {
            if record.is_some() && self.readers[index].advance() {
                heap.push(VcfRecordWithSource {
                    record: self.readers[index].current_record.clone(),
                    reader_index: index,
                });
            }
        }
        Ok(())
    }

    fn get_padding_base(&self, contig: &str, pos: i64) -> Vec<u8> {
        if let Some(genome_reader) = &self.genome_reader {
            genome_reader
                .fetch_seq(contig, pos as usize, pos as usize)
                .unwrap_or_else(|_| panic!("Failed to fetch sequence for chromosome {}", contig))
                .to_ascii_uppercase()
        } else {
            panic!("Genome reader is not available, but padding is required")
        }
    }

    fn add_padding_base(&mut self, sample_records: &mut [Option<Record>], contig: &str, pos: i64) {
        let padding_base = self.get_padding_base(contig, pos);

        for (record, reader) in sample_records.iter_mut().zip(&self.readers) {
            if reader.version.major >= Version::new(1, 0, 0).major {
                continue;
            }

            let Some(record) = record else { continue };

            let al_0 = record
                .format(b"AL")
                .integer()
                .expect("Error accessing FORMAT AL")[0]
                .iter()
                .min()
                .cloned()
                .unwrap();

            // Skip zero-length allele records
            if al_0 != 0 {
                let new_alleles: Vec<Vec<u8>> = record
                    .alleles()
                    .iter()
                    .map(|allele| {
                        let mut new_allele = Vec::with_capacity(1 + allele.len());
                        new_allele.extend_from_slice(&padding_base);
                        new_allele.extend_from_slice(allele);
                        new_allele
                    })
                    .collect();
                let new_alleles_refs: Vec<&[u8]> =
                    new_alleles.iter().map(|a| a.as_slice()).collect();
                record
                    .set_alleles(&new_alleles_refs)
                    .expect("Failed to set alleles")
            }
        }
    }

    fn merge_variant(
        &mut self,
        sample_records: &[Option<Record>],
        contig: &str,
        pos: i64,
    ) -> Result<()> {
        let template_record = self.get_template_record(sample_records);
        self.set_dummy_record_fields(template_record);

        let mut format_data = FormatData::new();
        // The outermost Vec represents different VCFs >  different samples within a VCF > genotype alleles of a sample
        let mut gt_vecs: Vec<Vec<Vec<GenotypeAllele>>> = Vec::new();
        let mut alleles: Vec<Vec<&[u8]>> = Vec::new();

        for (record_i, record) in sample_records.iter().enumerate() {
            if let Some(record) = record {
                self.process_record(record, &mut format_data, &mut gt_vecs, &mut alleles);
            } else {
                self.process_missing_record(record_i, &mut format_data, &mut gt_vecs, &mut alleles);
            }
        }

        self.merge_and_write_data(format_data, gt_vecs, alleles)
            .map_err(|e| format!("Failed to merge at {}:{}: {}.", contig, pos, e,))?;
        Ok(())
    }

    fn set_info_field(&mut self, record: &Record, field_name: &[u8], field_type: FieldType) {
        match field_type {
            FieldType::String => {
                let info_field = record.info(field_name).string().unwrap().unwrap();
                self.writer
                    .dummy_record
                    .push_info_string(field_name, &info_field)
                    .unwrap();
            }
            FieldType::Integer => {
                let info_field = record.info(field_name).integer().unwrap().unwrap();
                self.writer
                    .dummy_record
                    .push_info_integer(field_name, &info_field)
                    .unwrap();
            }
        }
    }

    fn get_template_record<'a>(&self, sample_records: &'a [Option<Record>]) -> &'a Record {
        let template_index = sample_records.iter().position(|r| r.is_some()).unwrap();
        sample_records[template_index].as_ref().unwrap()
    }

    fn set_dummy_record_fields(&mut self, template_record: &Record) {
        self.writer.dummy_record.set_rid(template_record.rid());
        self.writer.dummy_record.set_pos(template_record.pos());
        self.writer.dummy_record.set_qual(*MISSING_FLOAT);

        self.set_info_field(template_record, b"TRID", FieldType::String);
        self.set_info_field(template_record, b"END", FieldType::Integer);
        self.set_info_field(template_record, b"MOTIFS", FieldType::String);
        self.set_info_field(template_record, b"STRUC", FieldType::String);
    }

    fn process_record<'a>(
        &self,
        record: &'a Record,
        format_data: &mut FormatData,
        gt_vecs: &mut Vec<Vec<Vec<GenotypeAllele>>>,
        alleles: &mut Vec<Vec<&'a [u8]>>,
    ) {
        alleles.push(record.alleles().to_vec());
        self.process_genotypes(record, gt_vecs);
        self.process_format_fields(record, format_data);
    }

    fn process_missing_record(
        &self,
        record_i: usize,
        format_data: &mut FormatData,
        gt_vecs: &mut Vec<Vec<Vec<GenotypeAllele>>>,
        alleles: &mut Vec<Vec<&[u8]>>,
    ) {
        alleles.push(vec![]);
        let mut tmp_gt_vec = Vec::new();
        for _ in 0..self.readers[record_i].sample_n {
            tmp_gt_vec.push(vec![GenotypeAllele::UnphasedMissing]);
            PushMissingAndEnd::push_missing_and_end(&mut format_data.als);
            PushMissingAndEnd::push_missing_and_end(&mut format_data.sds);
            PushMissingAndEnd::push_missing_and_end(&mut format_data.aps);
            PushMissingAndEnd::push_missing_and_end(&mut format_data.ams);
            PushMissingAndEnd::push_missing_and_end(&mut format_data.allrs);
            PushMissingAndEnd::push_missing_and_end(&mut format_data.mcs);
            PushMissingAndEnd::push_missing_and_end(&mut format_data.mss);
        }
        gt_vecs.push(tmp_gt_vec);
    }

    fn process_genotypes(&self, record: &Record, gt_vecs: &mut Vec<Vec<Vec<GenotypeAllele>>>) {
        let gt_field = record.genotypes().unwrap();
        let mut file_gt_vecs: Vec<Vec<GenotypeAllele>> = Vec::new();
        for i in 0..record.sample_count() {
            let gt = gt_field.get(i as usize);
            file_gt_vecs.push(gt.iter().copied().collect());
        }
        gt_vecs.push(file_gt_vecs);
    }

    fn merge_and_write_data(
        &mut self,
        format_data: FormatData,
        gt_vecs: Vec<Vec<Vec<GenotypeAllele>>>,
        alleles: Vec<Vec<&[u8]>>,
    ) -> Result<()> {
        let (out_gts, out_alleles) = merge_exact(gt_vecs, alleles)?;
        self.writer
            .dummy_record
            .set_alleles(&out_alleles)
            .map_err(|e| e.to_string())?;

        let gts_new = self.flatten_genotypes(out_gts);
        self.write_format_fields(gts_new, format_data)?;

        self.writer
            .writer
            .write(&self.writer.dummy_record)
            .map_err(|e| e.to_string())?;
        self.writer.dummy_record.clear();
        Ok(())
    }

    fn flatten_genotypes(&self, out_gts: Vec<Vec<Vec<GenotypeAllele>>>) -> Vec<i32> {
        let mut gts_new = Vec::new();
        for file_gts in out_gts {
            for sample_gt in file_gts {
                let mut converted_sample_gt: Vec<i32> =
                    sample_gt.iter().map(|gt| i32::from(*gt)).collect();
                if converted_sample_gt.len() == 1 {
                    converted_sample_gt.push(VECTOR_END_INTEGER);
                }
                gts_new.extend(converted_sample_gt);
            }
        }
        gts_new
    }

    fn write_format_fields(&mut self, gts_new: Vec<i32>, format_data: FormatData) -> Result<()> {
        self.writer
            .dummy_record
            .push_format_integer(b"GT", &gts_new)
            .map_err(|e| e.to_string())?;
        self.writer
            .dummy_record
            .push_format_integer(b"AL", &format_data.als)
            .map_err(|e| e.to_string())?;
        self.writer
            .dummy_record
            .push_format_string(b"ALLR", &format_data.allrs)
            .map_err(|e| e.to_string())?;
        self.writer
            .dummy_record
            .push_format_integer(b"SD", &format_data.sds)
            .map_err(|e| e.to_string())?;
        self.writer
            .dummy_record
            .push_format_string(b"MC", &format_data.mcs)
            .map_err(|e| e.to_string())?;
        self.writer
            .dummy_record
            .push_format_string(b"MS", &format_data.mss)
            .map_err(|e| e.to_string())?;
        self.writer
            .dummy_record
            .push_format_float(b"AP", &format_data.aps)
            .map_err(|e| e.to_string())?;
        self.writer
            .dummy_record
            .push_format_float(b"AM", &format_data.ams)
            .map_err(|e| e.to_string())?;
        Ok(())
    }

    fn process_format_fields(&self, record: &Record, format_data: &mut FormatData) {
        self.process_integer_field(record, b"AL", &mut format_data.als);

        match record.format(b"ALLR").string() {
            Ok(_) => self.process_string_field(record, b"ALLR", &mut format_data.allrs),
            // Handle TRGT <=v0.3.4
            Err(_) => self.process_string_field(record, b"ALCI", &mut format_data.allrs),
        }

        self.process_integer_field(record, b"SD", &mut format_data.sds);
        self.process_string_field(record, b"MC", &mut format_data.mcs);
        self.process_string_field(record, b"MS", &mut format_data.mss);
        self.process_float_field(record, b"AP", &mut format_data.aps);
        self.process_am_field(record, &mut format_data.ams);
    }

    fn process_integer_field(&self, record: &Record, field_name: &[u8], values: &mut Vec<i32>) {
        let field = record.format(field_name).integer().unwrap_or_else(|_| {
            panic!(
                "Error accessing FORMAT {}",
                String::from_utf8_lossy(field_name)
            )
        });
        for sample_values in field.iter() {
            values.extend(sample_values.iter().copied());
            if sample_values.len() <= 1 {
                values.push(VECTOR_END_INTEGER);
            }
        }
    }

    fn process_float_field(&self, record: &Record, field_name: &[u8], values: &mut Vec<f32>) {
        let field = record.format(field_name).float().unwrap_or_else(|_| {
            panic!(
                "Error accessing FORMAT {}",
                String::from_utf8_lossy(field_name)
            )
        });
        for sample_values in field.iter() {
            values.extend(sample_values.iter().copied());
            if sample_values.len() <= 1 {
                values.push(*VECTOR_END_FLOAT);
            }
        }
    }

    fn process_string_field(&self, record: &Record, field_name: &[u8], values: &mut Vec<Vec<u8>>) {
        let field = record.format(field_name).string().unwrap_or_else(|_| {
            panic!(
                "Error accessing FORMAT {}",
                String::from_utf8_lossy(field_name)
            )
        });
        for sample_value in field.iter() {
            values.push(sample_value.to_vec());
        }
    }

    fn process_am_field(&self, record: &Record, ams: &mut Vec<f32>) {
        let am_field = record.format(b"AM");
        match am_field.float() {
            Ok(_) => self.process_float_field(record, b"AM", ams),
            // Handle TRGT <=v0.4.0
            Err(_) => {
                let int_field = record
                    .format(b"AM")
                    .integer()
                    .unwrap_or_else(|_| panic!("Error accessing FORMAT AM as an integer"));
                for sample_am in int_field.iter() {
                    let converted_am: Vec<f32> = sample_am
                        .iter()
                        .map(|&i| {
                            if i == i32::MIN {
                                *MISSING_FLOAT
                            } else {
                                i as f32 / 255.0
                            }
                        })
                        .collect();
                    ams.extend(converted_am);
                    if sample_am.len() <= 1 {
                        ams.push(*VECTOR_END_FLOAT);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf::{Read, Reader};
    use std::collections::BinaryHeap;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_vcf_record_wrapper_heap() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "##fileformat=VCFv4.2").unwrap();
        writeln!(temp_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        let reader = Reader::from_path(temp_file.path()).unwrap();

        let mut record1 = reader.empty_record();
        record1.set_rid(Some(1));
        record1.set_pos(100);

        let mut record2 = reader.empty_record();
        record2.set_rid(Some(1));
        record2.set_pos(2000);

        let mut record3 = reader.empty_record();
        record3.set_rid(Some(1));
        record3.set_pos(50);

        let mut record4 = reader.empty_record();
        record4.set_rid(Some(1));
        record4.set_pos(99);

        let mut heap = BinaryHeap::new();
        heap.push(VcfRecordWithSource {
            record: record1,
            reader_index: 0,
        });
        heap.push(VcfRecordWithSource {
            record: record4,
            reader_index: 3,
        });
        heap.push(VcfRecordWithSource {
            record: record2,
            reader_index: 1,
        });
        heap.push(VcfRecordWithSource {
            record: record3,
            reader_index: 2,
        });

        let r = heap.pop().unwrap();
        assert_eq!(r.record.rid(), Some(1));
        assert_eq!(r.record.pos(), 50);

        let r = heap.pop().unwrap();
        assert_eq!(r.record.rid(), Some(1));
        assert_eq!(r.record.pos(), 99);

        let r = heap.pop().unwrap();
        assert_eq!(r.record.rid(), Some(1));
        assert_eq!(r.record.pos(), 100);

        let r = heap.pop().unwrap();
        assert_eq!(r.record.rid(), Some(1));
        assert_eq!(r.record.pos(), 2000);
    }
}
