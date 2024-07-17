use super::{
    strategy::exact::merge_exact,
    vcf_reader::{VcfReader, VcfReaders},
    vcf_writer::VcfWriter,
};
use crate::{
    cli::MergeArgs,
    utils::{open_genome_reader, Result},
};
use once_cell::sync::Lazy;
use rust_htslib::{
    bcf::{self, header::HeaderView, record::GenotypeAllele, Record},
    faidx,
};
use semver::Version;
use std::{any::Any, cmp::Ordering, collections::BinaryHeap, env};

const MISSING_INTEGER: i32 = i32::MIN;
const VECTOR_END_INTEGER: i32 = i32::MIN + 1;
static MISSING_FLOAT: Lazy<f32> = Lazy::new(|| f32::from_bits(0x7F80_0001));
static VECTOR_END_FLOAT: Lazy<f32> = Lazy::new(|| f32::from_bits(0x7F80_0002));

fn _vec_to_comma_separated_string(vec: Vec<&[u8]>) -> String {
    vec.into_iter()
        .map(|slice| String::from_utf8_lossy(slice).to_string())
        .collect::<Vec<String>>()
        .join(",")
}

fn _header_to_string(header: &bcf::Header) -> String {
    unsafe {
        let header_ptr = header.inner;
        let mut header_len: i32 = 0;
        let header_cstr = rust_htslib::htslib::bcf_hdr_fmt_text(header_ptr, 0, &mut header_len);
        std::ffi::CStr::from_ptr(header_cstr)
            .to_string_lossy()
            .into_owned()
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
        self.record.rid() == other.record.rid() && self.record.pos() == other.record.pos()
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
        match self.record.rid().cmp(&other.record.rid()) {
            Ordering::Equal => self.record.pos().cmp(&other.record.pos()).reverse(),
            other => other.reverse(),
        }
    }
}

pub struct VcfProcessor {
    pub readers: Vec<VcfReader>,
    pub writer: VcfWriter,
    pub genome_reader: faidx::Reader, // TODO: Make optional? Only needed for <1.0
    // TODO: add args struct
    pub skip_n: usize,
    pub process_n: usize,
    pub needs_padding: bool,
}

impl VcfProcessor {
    pub fn new(args: &MergeArgs) -> Result<Self> {
        let vcf_readers = VcfReaders::new(args.vcfs.clone())?;
        if vcf_readers.readers.len() == 1 && !args.force_single {
            return Err("Expected two or more files to merge, got only one. Use --force-single to proceed anyway".into());
        }

        let genome_reader = open_genome_reader(&args.genome_path)?;
        let out_header = Self::create_output_header(&vcf_readers, args)?;
        let writer = VcfWriter::new(&out_header, &args.output_type, args.output.as_ref())?;

        let needs_padding = vcf_readers
            .readers
            .iter()
            .any(|reader| reader.version.major < Version::new(1, 0, 0).major);

        Ok(VcfProcessor {
            readers: vcf_readers.readers,
            writer,
            genome_reader,
            skip_n: args.skip_n.unwrap_or(0),
            process_n: args.process_n.unwrap_or(usize::MAX),
            needs_padding,
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

        let command_line = env::args().collect::<Vec<String>>().join(" ");
        let command_line = format!("##{}Command={}", env!("CARGO_PKG_NAME"), command_line);
        out_header.push_record(command_line.as_bytes());
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

    fn merge_variant(&mut self, sample_records: &[Option<Record>]) {
        let template_index = sample_records.iter().position(|r| r.is_some()).unwrap();
        let template_record = sample_records[template_index].as_ref().unwrap();

        self.writer.dummy_record.set_rid(template_record.rid());
        self.writer.dummy_record.set_pos(template_record.pos());
        self.writer.dummy_record.set_qual(template_record.qual());

        // TODO: Consolidate logic to allow for generic INFO fields: i32, f32, etc...
        self.set_info_field(template_record, b"TRID", FieldType::String);
        self.set_info_field(template_record, b"END", FieldType::Integer);
        self.set_info_field(template_record, b"MOTIFS", FieldType::String);
        self.set_info_field(template_record, b"STRUC", FieldType::String);

        // TODO: Clean this up
        // TODO: Consolidate logic to allow for generic FORMAT fields: i32, f32, etc...
        let mut als = Vec::new();
        let mut allrs = Vec::new();
        let mut sds = Vec::new();
        let mut mcs = Vec::new();
        let mut mss = Vec::new();
        let mut aps = Vec::new();
        let mut ams = Vec::new();
        let mut gt_vecs = Vec::new();
        let mut alleles = Vec::new();

        for record in sample_records.iter() {
            if let Some(record) = record {
                // TODO: Allow multiple Samples per record, at the moment we just take the first element
                alleles.push(record.alleles());

                // GT
                let gt_field = record.genotypes().unwrap();
                let gt = gt_field.get(0);
                gt_vecs.push(gt.iter().copied().collect());

                // TODO: Factor out redundancy
                let al_field = record
                    .format(b"AL")
                    .integer()
                    .expect("Error accessing FORMAT AL");
                als.extend(al_field[0].iter().copied());
                if al_field[0].len() == 1 {
                    als.push(VECTOR_END_INTEGER);
                }

                let allr = match record.format(b"ALLR").string() {
                    Ok(field) => field[0].to_vec(),
                    // Handle TRGT <=v0.3.4
                    Err(_) => {
                        let alci_field = record.format(b"ALCI").string().unwrap();
                        alci_field[0].to_vec()
                    }
                };
                allrs.push(allr);

                let sd_field = record
                    .format(b"SD")
                    .integer()
                    .expect("Error accessing FORMAT SD");
                sds.extend(sd_field[0].iter().copied());
                if sd_field[0].len() == 1 {
                    sds.push(VECTOR_END_INTEGER);
                }

                let mc_field = record
                    .format(b"MC")
                    .string()
                    .expect("Error acessing FORMAT MC");
                let mc = mc_field[0].to_vec();
                mcs.push(mc);

                let ms_field = record
                    .format(b"MS")
                    .string()
                    .expect("Error acessing FORMAT MS");
                let ms = ms_field[0].to_vec();
                mss.push(ms);

                let ap_field = record
                    .format(b"AP")
                    .float()
                    .expect("Error accessing FORMAT AP");
                aps.extend(ap_field[0].iter().copied());
                if ap_field[0].len() == 1 {
                    aps.push(*VECTOR_END_FLOAT);
                }

                let am_field = match record.format(b"AM").float() {
                    Ok(field) => field[0].to_vec(),
                    // Handle TRGT <=v0.4.0
                    Err(_) => {
                        let int_field = record
                            .format(b"AM")
                            .integer()
                            .expect("Error accessing FORMAT AM as an integer");
                        int_field[0]
                            .iter()
                            .map(|&i| {
                                // Account for missing values
                                if i == i32::MIN {
                                    *MISSING_FLOAT
                                } else {
                                    i as f32 / 255.0
                                }
                            })
                            .collect::<Vec<_>>()
                    }
                };
                ams.extend(am_field.iter().copied());
                if am_field.len() == 1 {
                    ams.push(*VECTOR_END_FLOAT);
                }
            } else {
                gt_vecs.push(vec![GenotypeAllele::UnphasedMissing]);
                alleles.push(vec![]);

                PushMissingAndEnd::push_missing_and_end(&mut als);
                PushMissingAndEnd::push_missing_and_end(&mut sds);
                PushMissingAndEnd::push_missing_and_end(&mut aps);
                PushMissingAndEnd::push_missing_and_end(&mut ams);
                PushMissingAndEnd::push_missing_and_end(&mut allrs);
                PushMissingAndEnd::push_missing_and_end(&mut mcs);
                PushMissingAndEnd::push_missing_and_end(&mut mss);
            }
        }

        // Merge alleles and genotypes
        let (out_gts, out_alleles) = merge_exact(gt_vecs, alleles);
        self.writer.dummy_record.set_alleles(&out_alleles).unwrap();

        // Flatten to a 1D 2D representation using
        let mut gts_new: Vec<i32> = Vec::new();
        for sample_gt in out_gts {
            let mut converted_sample_gt: Vec<i32> =
                sample_gt.iter().map(|gt| i32::from(*gt)).collect();
            if converted_sample_gt.len() == 1 {
                converted_sample_gt.push(VECTOR_END_INTEGER);
            }
            gts_new.extend(converted_sample_gt);
        }
        //

        self.writer
            .dummy_record
            .push_format_integer(b"GT", &gts_new)
            .unwrap();

        self.writer
            .dummy_record
            .push_format_integer(b"AL", &als)
            .unwrap();

        self.writer
            .dummy_record
            .push_format_string(b"ALLR", &allrs)
            .unwrap();

        self.writer
            .dummy_record
            .push_format_integer(b"SD", &sds)
            .unwrap();

        self.writer
            .dummy_record
            .push_format_string(b"MC", &mcs)
            .unwrap();

        self.writer
            .dummy_record
            .push_format_string(b"MS", &mss)
            .unwrap();

        self.writer
            .dummy_record
            .push_format_float(b"AP", &aps)
            .unwrap();

        self.writer
            .dummy_record
            .push_format_float(b"AM", &ams)
            .unwrap();

        self.writer.writer.write(&self.writer.dummy_record).unwrap();

        self.writer.dummy_record.clear();
    }

    fn init_heap(&mut self) -> BinaryHeap<VcfRecordWithSource> {
        let mut heap = BinaryHeap::new();
        for (index, reader) in self.readers.iter_mut().enumerate() {
            if reader.advance() {
                heap.push(VcfRecordWithSource {
                    record: reader.current_record.clone(),
                    reader_index: index,
                });
            }
        }
        heap
    }

    fn update_heap(
        &mut self,
        heap: &mut BinaryHeap<VcfRecordWithSource>,
        sample_records: &[Option<Record>],
    ) {
        for (index, record) in sample_records.iter().enumerate() {
            if record.is_some() && self.readers[index].advance() {
                heap.push(VcfRecordWithSource {
                    record: self.readers[index].current_record.clone(),
                    reader_index: index,
                });
            }
        }
    }

    pub fn merge_variants(&mut self) {
        let mut n = 0;
        let mut n_processed = 0;

        let mut sample_records = vec![None; self.readers.len()];
        let mut heap = self.init_heap();
        while let Some(min_element) = heap.pop() {
            let min_rid = min_element.record.rid().unwrap();
            let min_pos = min_element.record.pos();
            sample_records[min_element.reader_index] = Some(min_element.record);

            while let Some(peek_next_element) = heap.peek() {
                if peek_next_element.record.rid().unwrap() == min_rid
                    && peek_next_element.record.pos() == min_pos
                {
                    let next_element = heap.pop().unwrap();
                    sample_records[next_element.reader_index] = Some(next_element.record);
                } else {
                    break;
                }
            }

            if n >= self.skip_n {
                log::info!("Processing: {}:{}", min_rid, min_pos);
                if self.needs_padding {
                    let padding_base = self.get_padding_base(
                        min_rid,
                        min_pos,
                        &self.readers[min_element.reader_index].header,
                    );
                    self.add_padding_base(&mut sample_records, padding_base);
                }

                self.merge_variant(&sample_records);
                n_processed += 1;
                if n_processed >= self.process_n {
                    break;
                }
            }
            n += 1;

            self.update_heap(&mut heap, &sample_records);
            sample_records.fill(None);
        }
    }

    fn add_padding_base(&mut self, sample_records: &mut [Option<Record>], padding_base: Vec<u8>) {
        for (index, record) in sample_records.iter_mut().enumerate() {
            if self.readers[index].version.major < Version::new(1, 0, 0).major {
                if let Some(record) = record {
                    let al_0 = record
                        .format(b"AL")
                        .integer()
                        .expect("Error accessing FORMAT AL")[0]
                        .iter()
                        .min()
                        .cloned()
                        .unwrap();
                    // Zero-length allele records do not need to be updated
                    if al_0 != 0 {
                        let new_alleles: Vec<Vec<u8>> = record
                            .alleles()
                            .iter()
                            .map(|allele| {
                                let mut new_allele = padding_base.to_vec();
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
        }
    }

    fn get_padding_base(&self, rid: u32, pos: i64, header: &HeaderView) -> Vec<u8> {
        let chrom = header.rid2name(rid).unwrap();
        let chrom_str = std::str::from_utf8(chrom).expect("Invalid UTF-8 sequence");
        self.genome_reader
            .fetch_seq(chrom_str, pos as usize, pos as usize)
            .map(|seq| seq.to_vec())
            .ok()
            .unwrap()
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
        record4.set_rid(Some(10));
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
        assert_eq!(r.record.pos(), 100);

        let r = heap.pop().unwrap();
        assert_eq!(r.record.rid(), Some(1));
        assert_eq!(r.record.pos(), 2000);

        let r = heap.pop().unwrap();
        assert_eq!(r.record.rid(), Some(10));
        assert_eq!(r.record.pos(), 99);
    }
}
