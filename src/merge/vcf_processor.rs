use super::{
    field_descriptors::{FORMAT_FIELDS, INFO_FIELDS},
    field_registry::{
        build_predefined_field_ids, FieldCategory, FieldDataRegistry, FieldDescriptor,
        FieldIdCache, VcfType, MISSING_FLOAT, VECTOR_END_INTEGER,
    },
    strategy::exact::merge_exact,
    tpool::HtsThreadPool,
    vcf_reader::{VcfReadMode, VcfReaders},
    vcf_writer::VcfWriter,
};
use crate::{
    cli::MergeArgs,
    utils::{format_number_with_commas, open_genome_reader, Result},
};
use rust_htslib::{
    bcf::{self, record::GenotypeAllele, HeaderRecord, Record},
    faidx,
};
use std::{
    cmp::Ordering,
    collections::{BTreeMap, BinaryHeap, HashMap, HashSet},
    env, mem,
    path::PathBuf,
    sync::Arc,
};

struct CollectedVariantData<'a> {
    /// Genotypes, organized as: VCF -> sample -> genotype alleles
    genotypes: Vec<Vec<Vec<GenotypeAllele>>>,
    /// Alleles from each VCF
    alleles: Vec<Vec<&'a [u8]>>,
}

type TridGroups = BTreeMap<Vec<u8>, Vec<(usize, Record)>>;

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

#[derive(Debug)]
struct StreamRecordWithSource {
    record: bcf::Record,
    reader_index: usize,
    contig_idx: usize,
}

impl PartialEq for StreamRecordWithSource {
    fn eq(&self, other: &Self) -> bool {
        self.contig_idx == other.contig_idx && self.record.pos() == other.record.pos()
    }
}

impl Eq for StreamRecordWithSource {}

impl PartialOrd for StreamRecordWithSource {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for StreamRecordWithSource {
    fn cmp(&self, other: &Self) -> Ordering {
        self.contig_idx
            .cmp(&other.contig_idx)
            .then_with(|| self.record.pos().cmp(&other.record.pos()))
            .reverse()
    }
}

#[derive(Default)]
struct MergeProgress {
    seen: usize,
    processed: usize,
    failed: usize,
}

pub struct VcfProcessor {
    // writer must drop BEFORE vcf_readers to properly flush BGZF buffer
    writer: VcfWriter,
    vcf_readers: VcfReaders,
    genome_reader: Option<faidx::Reader>,
    contig_order: Vec<String>,
    contig_to_idx: HashMap<String, usize>,
    read_mode: VcfReadMode,
    skip_n: usize,
    process_n: usize,
    needs_padding: bool,
    quit_on_error: bool,
    active_format_fields: Vec<FieldDescriptor<'static>>,
    active_info_fields: Vec<FieldDescriptor<'static>>,
    field_registry: FieldDataRegistry,
    field_id_caches: Vec<FieldIdCache>,
}

impl VcfProcessor {
    pub fn new(args: &MergeArgs, vcfs: Vec<PathBuf>) -> Result<Self> {
        let read_mode = args.no_index;

        let thread_pool = if args.threads > 1 {
            log::debug!("Building shared thread pool with {} threads", args.threads);
            Some(Arc::new(HtsThreadPool::new(args.threads as i32)?))
        } else {
            None
        };

        let vcf_readers = VcfReaders::new(vcfs, read_mode, thread_pool.clone())?;
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

        let contig_to_idx: HashMap<String, usize> = contig_order
            .iter()
            .enumerate()
            .map(|(idx, contig)| (contig.clone(), idx))
            .collect();

        let needs_padding = vcf_readers
            .readers
            .iter()
            .any(|reader| reader.needs_pos_adjustment);

        let genome_reader = if needs_padding {
            Some(open_genome_reader(args.genome_src.as_ref().ok_or(
                "A reference genome is required for merging pre v1.0 TRGT VCFs, provide as --genome ref.fa"
            )?)?)
        } else {
            None
        };

        let out_header = Self::create_output_header(&vcf_readers, args)?;
        let writer = VcfWriter::new(
            &out_header,
            args.output_type.as_ref(),
            args.output.as_ref(),
            thread_pool,
        )?;

        if needs_padding {
            log::debug!("At least one VCF file is pre-1.0 and needs base padding!");
        }

        let num_readers = vcf_readers.readers.len();
        let mut processor = VcfProcessor {
            writer,
            vcf_readers,
            genome_reader,
            contig_order,
            contig_to_idx,
            read_mode,
            skip_n: args.skip_n.unwrap_or(0),
            process_n: args.process_n.unwrap_or(usize::MAX),
            needs_padding,
            quit_on_error: args.quit_on_error,
            active_format_fields: Vec::new(),
            active_info_fields: Vec::new(),
            field_registry: FieldDataRegistry::from_descriptors(&[]),
            field_id_caches: (0..num_readers).map(|_| FieldIdCache::new()).collect(),
        };

        processor.discover_fields()?;
        processor.field_registry =
            FieldDataRegistry::from_descriptors(&processor.active_format_fields);

        Ok(processor)
    }

    fn create_output_header(vcf_readers: &VcfReaders, args: &MergeArgs) -> Result<bcf::Header> {
        let mut out_header = bcf::Header::new();
        vcf_readers.merge_headers(&mut out_header)?;

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
            crate::cli::FULL_VERSION
        );
        out_header.push_record(version_line.as_bytes());

        let command_line = format!(
            "##{}Command={}",
            env!("CARGO_PKG_NAME"),
            env::args().collect::<Vec<String>>().join(" ")
        );
        out_header.push_record(command_line.as_bytes());
    }

    fn discover_fields(&mut self) -> Result<()> {
        self.active_info_fields =
            self.discover_predefined_fields(INFO_FIELDS, FieldCategory::Info)?;
        self.active_info_fields
            .extend(self.discover_custom_fields(FieldCategory::Info, INFO_FIELDS));

        self.active_format_fields =
            self.discover_predefined_fields(FORMAT_FIELDS, FieldCategory::Format)?;
        self.active_format_fields
            .extend(self.discover_custom_fields(FieldCategory::Format, FORMAT_FIELDS));

        self.log_discovery_summary();
        Ok(())
    }

    fn discover_predefined_fields(
        &self,
        descriptors: &'static [FieldDescriptor<'static>],
        category: FieldCategory,
    ) -> Result<Vec<FieldDescriptor<'static>>> {
        let mut active_ids: HashSet<&'static [u8]> = HashSet::new();
        let category_name = match category {
            FieldCategory::Format => "FORMAT",
            FieldCategory::Info => "INFO",
        };

        for descriptor in descriptors {
            // Skip special handling fields for FORMAT category
            if category == FieldCategory::Format && descriptor.special_handling {
                continue;
            }

            let present_count = self
                .vcf_readers
                .readers
                .iter()
                .filter(|reader| Self::field_exists_in_header(&reader.header, descriptor))
                .count();

            if present_count > 0 {
                let id = descriptor
                    .borrowed_id()
                    .expect("Predefined descriptors must have borrowed IDs");
                active_ids.insert(id);
            }

            if descriptor.required && present_count != self.vcf_readers.len() {
                return Err(format!(
                    "Required {} field {} is missing from {}/{} input VCF(s)",
                    category_name,
                    String::from_utf8_lossy(descriptor.id.as_ref()),
                    self.vcf_readers.len() - present_count,
                    self.vcf_readers.len()
                ));
            }
        }

        Ok(descriptors
            .iter()
            .filter(|d| {
                d.borrowed_id()
                    .map(|id| active_ids.contains(&id))
                    .unwrap_or(false)
            })
            .cloned()
            .collect())
    }

    fn validate_field_consistency(
        field_id: &[u8],
        descriptors: &[FieldDescriptor<'static>],
        category: FieldCategory,
    ) -> Option<FieldDescriptor<'static>> {
        let first = &descriptors[0];
        let all_consistent = descriptors
            .iter()
            .all(|d| d.vcf_type == first.vcf_type && d.number == first.number);

        if !all_consistent {
            let category_name = match category {
                FieldCategory::Format => "FORMAT",
                FieldCategory::Info => "INFO",
            };
            log::warn!(
                "Custom {} field {} has inconsistent definitions across VCFs. Skipping.",
                category_name,
                String::from_utf8_lossy(field_id)
            );
            return None;
        }

        log::debug!(
            "Discovered custom {} field: {}",
            match category {
                FieldCategory::Format => "FORMAT",
                FieldCategory::Info => "INFO",
            },
            String::from_utf8_lossy(field_id)
        );
        Some(first.clone())
    }

    fn discover_custom_fields(
        &self,
        category: FieldCategory,
        predefined_descriptors: &'static [FieldDescriptor<'static>],
    ) -> Vec<FieldDescriptor<'static>> {
        let predefined_ids = build_predefined_field_ids(predefined_descriptors);
        let mut custom_fields: BTreeMap<Vec<u8>, Vec<FieldDescriptor<'static>>> = BTreeMap::new();

        for reader in &self.vcf_readers.readers {
            for record in reader.header.header_records() {
                let values = match (record, category) {
                    (HeaderRecord::Info { values, .. }, FieldCategory::Info) => values,
                    (HeaderRecord::Format { values, .. }, FieldCategory::Format) => values,
                    _ => continue,
                };

                let Some(id) = values.get("ID") else { continue };
                let Some(type_str) = values.get("Type") else {
                    continue;
                };
                let Some(number_str) = values.get("Number") else {
                    continue;
                };
                let description = values.get("Description").map(|s| s.as_str()).unwrap_or("");

                if let Some(descriptor) = FieldDescriptor::from_header_values(
                    id,
                    type_str,
                    number_str,
                    description,
                    category,
                ) {
                    if predefined_ids.contains(descriptor.id.as_ref()) {
                        continue;
                    }
                    custom_fields
                        .entry(descriptor.id.as_ref().to_vec())
                        .or_default()
                        .push(descriptor);
                }
            }
        }

        custom_fields
            .into_iter()
            .filter_map(|(field_id, descriptors)| {
                Self::validate_field_consistency(&field_id, &descriptors, category)
            })
            .collect()
    }

    fn log_discovery_summary(&self) {
        let predefined_format_count = FORMAT_FIELDS.iter().filter(|d| !d.special_handling).count();
        let predefined_info_count = INFO_FIELDS.len();

        log::debug!(
            "Discovered {} active FORMAT fields ({} predefined, {} custom) and {} active INFO fields ({} predefined, {} custom)",
            self.active_format_fields.len(),
            predefined_format_count,
            self.active_format_fields.len().saturating_sub(predefined_format_count),
            self.active_info_fields.len(),
            predefined_info_count,
            self.active_info_fields.len().saturating_sub(predefined_info_count)
        );
    }

    fn field_exists_in_header(
        header: &bcf::header::HeaderView,
        descriptor: &FieldDescriptor<'static>,
    ) -> bool {
        let category_key = match descriptor.category {
            FieldCategory::Format => "FORMAT",
            FieldCategory::Info => "INFO",
        };

        for record in header.header_records() {
            match record {
                HeaderRecord::Format { values, .. } if category_key == "FORMAT" => {
                    if let Some(id) = values.get("ID") {
                        if id.as_bytes() == descriptor.id.as_ref() {
                            return true;
                        }
                        for &alias in descriptor.aliases {
                            if id.as_bytes() == alias {
                                return true;
                            }
                        }
                    }
                }
                HeaderRecord::Info { values, .. } if category_key == "INFO" => {
                    if let Some(id) = values.get("ID") {
                        if id.as_bytes() == descriptor.id.as_ref() {
                            return true;
                        }
                        for &alias in descriptor.aliases {
                            if id.as_bytes() == alias {
                                return true;
                            }
                        }
                    }
                }
                _ => {}
            }
        }
        false
    }

    pub fn merge_variants(&mut self) -> Result<()> {
        log::debug!("Merging variants in {} mode", self.read_mode);
        match self.read_mode {
            VcfReadMode::Indexed => self.merge_variants_indexed(),
            VcfReadMode::Streaming => self.merge_variants_streaming(),
        }
    }

    fn merge_variants_indexed(&mut self) -> Result<()> {
        let mut progress = MergeProgress::default();

        let mut sample_records = vec![None; self.vcf_readers.len()];
        let mut per_reader_records: Vec<Vec<Record>> = vec![Vec::new(); self.vcf_readers.len()];
        let mut future_entries: Vec<VcfRecordWithSource> = Vec::new();

        for contig in self.contig_order.clone() {
            let mut heap = self.init_heap(&contig)?;

            while let Some(min_element) = heap.pop() {
                future_entries.clear();

                let min_pos = min_element.record.pos();
                per_reader_records[min_element.reader_index].push(min_element.record);

                while let Some(peek_next_element) = heap.peek() {
                    if peek_next_element.record.pos() == min_pos {
                        let next_element = heap.pop().unwrap();
                        per_reader_records[next_element.reader_index].push(next_element.record);
                    } else {
                        break;
                    }
                }

                for (reader_idx, records) in per_reader_records.iter_mut().enumerate() {
                    if records.is_empty() {
                        continue;
                    }
                    loop {
                        if !self.vcf_readers.readers[reader_idx].advance() {
                            break;
                        }
                        let record = self.vcf_readers.readers[reader_idx].take_current_record();
                        if record.pos() == min_pos {
                            records.push(record);
                        } else {
                            future_entries.push(VcfRecordWithSource {
                                record,
                                reader_index: reader_idx,
                            });
                            break;
                        }
                    }
                }

                for entry in future_entries.drain(..) {
                    heap.push(entry);
                }

                let grouped =
                    self.group_records_by_trid(&mut per_reader_records, &contig, min_pos)?;
                if self.process_grouped_records(
                    grouped,
                    &contig,
                    min_pos,
                    &mut sample_records,
                    &mut progress,
                )? {
                    return Ok(());
                }
            }
        }
        self.log_merge_summary(&progress);
        Ok(())
    }

    fn merge_variants_streaming(&mut self) -> Result<()> {
        let mut progress = MergeProgress::default();

        let mut sample_records = vec![None; self.vcf_readers.len()];
        let mut per_reader_records: Vec<Vec<Record>> = vec![Vec::new(); self.vcf_readers.len()];
        let mut future_entries: Vec<StreamRecordWithSource> = Vec::new();

        // Track contig and position per reader to enforce streaming order without random access
        let mut last_rid_per_reader: Vec<Option<u32>> = vec![None; self.vcf_readers.len()];
        let mut last_pos_per_reader: Vec<Option<i64>> = vec![None; self.vcf_readers.len()];
        let mut heap =
            self.init_streaming_heap(&mut last_rid_per_reader, &mut last_pos_per_reader)?;

        while let Some(min_element) = heap.pop() {
            future_entries.clear();

            let contig = self
                .contig_order
                .get(min_element.contig_idx)
                .ok_or_else(|| format!("Invalid contig index {}", min_element.contig_idx))?
                .clone();
            let min_contig_idx = min_element.contig_idx;
            let min_pos = min_element.record.pos();
            per_reader_records[min_element.reader_index].push(min_element.record);

            while let Some(peek_next_element) = heap.peek() {
                if peek_next_element.contig_idx == min_element.contig_idx
                    && peek_next_element.record.pos() == min_pos
                {
                    let next_element = heap.pop().unwrap();
                    per_reader_records[next_element.reader_index].push(next_element.record);
                } else {
                    break;
                }
            }

            for (reader_idx, records) in per_reader_records.iter_mut().enumerate() {
                if records.is_empty() {
                    continue;
                }
                loop {
                    match self.read_next_stream_record(
                        reader_idx,
                        &mut last_rid_per_reader,
                        &mut last_pos_per_reader,
                    )? {
                        Some(next_entry)
                            if next_entry.contig_idx == min_contig_idx
                                && next_entry.record.pos() == min_pos =>
                        {
                            records.push(next_entry.record);
                        }
                        Some(next_entry) => {
                            future_entries.push(next_entry);
                            break;
                        }
                        None => break,
                    }
                }
            }

            for entry in future_entries.drain(..) {
                heap.push(entry);
            }

            let grouped = self.group_records_by_trid(&mut per_reader_records, &contig, min_pos)?;
            if self.process_grouped_records(
                grouped,
                &contig,
                min_pos,
                &mut sample_records,
                &mut progress,
            )? {
                return Ok(());
            }
        }

        self.log_merge_summary(&progress);
        Ok(())
    }

    fn process_grouped_records(
        &mut self,
        grouped: TridGroups,
        contig: &str,
        pos: i64,
        sample_records: &mut [Option<Record>],
        progress: &mut MergeProgress,
    ) -> Result<bool> {
        for (trid, group_records) in grouped {
            for slot in sample_records.iter_mut() {
                *slot = None;
            }

            self.load_group_into_sample_records(sample_records, group_records, contig, pos, &trid)?;

            let should_process = progress.seen >= self.skip_n;
            let merge_result = if should_process {
                log::trace!(
                    "Processing: {}:{} (TRID={})",
                    contig,
                    pos,
                    String::from_utf8_lossy(&trid)
                );
                if self.needs_padding {
                    self.add_padding_base(sample_records, contig, pos);
                }
                self.merge_variant(sample_records, contig, pos)
            } else {
                Ok(())
            };

            self.return_records(sample_records);

            if should_process {
                match merge_result {
                    Ok(()) => {
                        progress.processed += 1;
                        if progress.processed >= self.process_n {
                            return Ok(true);
                        }
                    }
                    Err(e) => {
                        if self.quit_on_error {
                            return Err(e);
                        }
                        progress.failed += 1;
                        log::warn!("{}. Skipping...", e);
                    }
                }
            } else if let Err(e) = merge_result {
                if self.quit_on_error {
                    return Err(e);
                }
            }

            progress.seen += 1;
        }
        Ok(false)
    }

    fn log_merge_summary(&self, progress: &MergeProgress) {
        let mut log_message = format!(
            "Successfully merged {} TR site(s).",
            format_number_with_commas(progress.processed)
        );
        if progress.failed > 0 {
            log_message.push_str(&format!(
                " Failed to merge {} TR site(s)!",
                format_number_with_commas(progress.failed)
            ));
        }
        log::info!("{}", log_message);
    }

    fn init_heap(&mut self, contig: &str) -> Result<BinaryHeap<VcfRecordWithSource>> {
        let mut heap = BinaryHeap::new();
        for (index, reader) in self.vcf_readers.readers.iter_mut().enumerate() {
            let rid = match reader.header.name2rid(contig.as_bytes()) {
                Ok(id) => id,
                Err(_) => continue,
            };

            if let Some(indexed_reader) = reader.indexed_reader_mut() {
                if indexed_reader.fetch(rid, 0, None).is_err() {
                    continue;
                }

                if reader.advance() {
                    heap.push(VcfRecordWithSource {
                        record: reader.take_current_record(),
                        reader_index: index,
                    });
                }
            }
        }
        Ok(heap)
    }

    fn init_streaming_heap(
        &mut self,
        last_rid_per_reader: &mut [Option<u32>],
        last_pos_per_reader: &mut [Option<i64>],
    ) -> Result<BinaryHeap<StreamRecordWithSource>> {
        let mut heap = BinaryHeap::new();
        for reader_idx in 0..self.vcf_readers.readers.len() {
            if let Some(entry) =
                self.read_next_stream_record(reader_idx, last_rid_per_reader, last_pos_per_reader)?
            {
                heap.push(entry);
            }
        }
        Ok(heap)
    }

    fn read_next_stream_record(
        &mut self,
        reader_idx: usize,
        last_rid_per_reader: &mut [Option<u32>],
        last_pos_per_reader: &mut [Option<i64>],
    ) -> Result<Option<StreamRecordWithSource>> {
        loop {
            if !self.vcf_readers.readers[reader_idx].advance() {
                return Ok(None);
            }
            let record = self.vcf_readers.readers[reader_idx].take_current_record();
            let pos = record.pos();
            let rid = record.rid().ok_or_else(|| {
                format!(
                    "Record missing RID in {}",
                    self.vcf_readers.readers[reader_idx].file_path
                )
            })?;

            if let Some(prev_rid) = last_rid_per_reader[reader_idx] {
                match rid.cmp(&prev_rid) {
                    Ordering::Less => {
                        let header = &self.vcf_readers.readers[reader_idx].header;
                        let current_name = header
                            .rid2name(rid)
                            .map(|n| String::from_utf8_lossy(n).into_owned())
                            .unwrap_or_else(|_| "<unknown>".to_string());
                        let previous_name = header
                            .rid2name(prev_rid)
                            .map(|n| String::from_utf8_lossy(n).into_owned())
                            .unwrap_or_else(|_| "<unknown>".to_string());
                        return Err(format!(
                            "Input VCF {} is not sorted correctly for streaming mode: contig {} appears after {}.\n\
                        For streaming merge (--no-index), all VCFs must have identical contig order and be sorted.\n\
                        Either: (1) re-sort this VCF to match others, or (2) use indexed mode (remove --no-index)",
                           self.vcf_readers.readers[reader_idx].file_path, current_name, previous_name
                        ));
                    }
                    Ordering::Equal => {
                        if let Some(prev_pos) = last_pos_per_reader[reader_idx] {
                            if pos < prev_pos {
                                let header = &self.vcf_readers.readers[reader_idx].header;
                                let contig_name = header
                                    .rid2name(rid)
                                    .map(|n| String::from_utf8_lossy(n).into_owned())
                                    .unwrap_or_else(|_| "<unknown>".to_string());
                                return Err(format!(
                                    "Input VCF {} is not sorted correctly within contig {}: position {} appears after {}.\n\
                        For streaming merge (--no-index), each contig must be position-sorted.\n\
                        Either: (1) re-sort this VCF, or (2) use indexed mode (remove --no-index)",
                                    self.vcf_readers.readers[reader_idx].file_path,
                                    contig_name,
                                    pos + 1,
                                    prev_pos + 1
                                ));
                            }
                        }
                    }
                    Ordering::Greater => {
                        last_pos_per_reader[reader_idx] = None;
                    }
                }
            }
            last_rid_per_reader[reader_idx] = Some(rid);
            last_pos_per_reader[reader_idx] = Some(pos);

            let contig_bytes = self.vcf_readers.readers[reader_idx]
                .header
                .rid2name(rid)
                .map_err(|e| e.to_string())?;
            let contig_name = String::from_utf8_lossy(contig_bytes).into_owned();

            let Some(&contig_idx) = self.contig_to_idx.get(&contig_name) else {
                self.vcf_readers.readers[reader_idx].return_record(record);
                continue;
            };

            return Ok(Some(StreamRecordWithSource {
                record,
                reader_index: reader_idx,
                contig_idx,
            }));
        }
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

        for (record, reader) in sample_records.iter_mut().zip(&self.vcf_readers.readers) {
            if !reader.needs_pos_adjustment {
                continue;
            }

            let Some(record) = record else { continue };

            let al_fmt = record
                .format(b"AL")
                .integer()
                .expect("Error accessing FORMAT AL");
            let al_0 = *al_fmt[0].iter().min().unwrap();

            // Zero-length allele records already have a padding base, so skip them
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
        self.field_registry.clear();
        let collected = self
            .collect_variant_data(sample_records)
            .map_err(|e| format!("Failed to collect data at {}:{}: {}", contig, pos + 1, e))?;

        self.merge_and_write(collected, sample_records)
            .map_err(|e| format!("Failed to merge/write at {}:{}: {}", contig, pos + 1, e))?;

        Ok(())
    }

    fn group_records_by_trid(
        &self,
        per_reader_records: &mut [Vec<Record>],
        contig: &str,
        pos: i64,
    ) -> Result<TridGroups> {
        let mut grouped: TridGroups = BTreeMap::new();
        for (reader_idx, records) in per_reader_records.iter_mut().enumerate() {
            for record in records.drain(..) {
                let trid = Self::extract_trid(&record, contig, pos)?;
                grouped.entry(trid).or_default().push((reader_idx, record));
            }
        }
        Ok(grouped)
    }

    fn load_group_into_sample_records(
        &self,
        sample_records: &mut [Option<Record>],
        mut group_records: Vec<(usize, Record)>,
        contig: &str,
        pos: i64,
        trid: &[u8],
    ) -> Result<()> {
        for (reader_idx, record) in group_records.drain(..) {
            if reader_idx >= sample_records.len() {
                return Err(format!(
                    "Reader index {} is out of bounds when processing {}:{}",
                    reader_idx,
                    contig,
                    pos + 1
                ));
            }
            if sample_records[reader_idx].is_some() {
                let reader_name = self
                    .vcf_readers
                    .readers
                    .get(reader_idx)
                    .map(|reader| reader.file_path.as_str())
                    .unwrap_or("<unknown VCF>");
                return Err(format!(
                    "Multiple records detected for {}:{} (TRID={}) in reader: {}",
                    contig,
                    pos + 1,
                    String::from_utf8_lossy(trid),
                    reader_name
                ));
            }
            sample_records[reader_idx] = Some(record);
        }
        Ok(())
    }

    fn return_records(&mut self, sample_records: &mut [Option<Record>]) {
        for (reader_idx, record_slot) in sample_records.iter_mut().enumerate() {
            if let Some(record) = record_slot.take() {
                self.vcf_readers.readers[reader_idx].return_record(record);
            }
        }
    }

    fn extract_trid(record: &Record, contig: &str, pos: i64) -> Result<Vec<u8>> {
        let values = record
            .info(b"TRID")
            .string()
            .map_err(|e| format!("Failed to read TRID at {}:{}: {}", contig, pos + 1, e))?
            .ok_or_else(|| format!("Missing TRID at {}:{}", contig, pos + 1))?;
        if values.len() != 1 {
            return Err(format!(
                "Expected exactly one TRID value at {}:{}, found {}",
                contig,
                pos + 1,
                values.len()
            ));
        }
        Ok(values[0].to_vec())
    }

    fn collect_variant_data<'a>(
        &mut self,
        sample_records: &'a [Option<Record>],
    ) -> Result<CollectedVariantData<'a>> {
        let n_vcfs = sample_records.len();
        let mut genotypes: Vec<Vec<Vec<GenotypeAllele>>> = Vec::with_capacity(n_vcfs);
        let mut alleles: Vec<Vec<&[u8]>> = Vec::with_capacity(n_vcfs);

        for (record_i, record_opt) in sample_records.iter().enumerate() {
            if let Some(record) = record_opt {
                let file_path = self
                    .vcf_readers
                    .readers
                    .get(record_i)
                    .map(|r| r.file_path.as_str())
                    .unwrap_or("<unknown source>");

                Self::check_required_info_fields(
                    record,
                    file_path,
                    &self.active_info_fields,
                    &mut self.field_id_caches[record_i],
                )?;

                alleles.push(record.alleles());
                if let Err(e) = Self::process_genotypes(record, &mut genotypes) {
                    return Err(format!(
                        "Error processing GT field for file {}: {}",
                        file_path, e
                    ));
                }
                self.field_registry
                    .append_record(
                        record,
                        &self.active_format_fields,
                        &mut self.field_id_caches[record_i],
                    )
                    .map_err(|e| format!("{} has invalid FORMAT field. {}", file_path, e))?;
            } else {
                // Add missing alleles, genotypes, and FORMAT field data for VCF file that have no record at this position
                alleles.push(vec![]);
                let sample_count = self.vcf_readers.readers[record_i].sample_n as usize;
                let file_genotypes = (0..sample_count)
                    .map(|_| vec![GenotypeAllele::UnphasedMissing])
                    .collect::<Vec<_>>();
                genotypes.push(file_genotypes);
                self.field_registry
                    .push_missing_for_all_fields(sample_count);
            }
        }

        Ok(CollectedVariantData { genotypes, alleles })
    }

    fn check_required_info_fields(
        record: &Record,
        reader_file_path: &str,
        info_fields: &[FieldDescriptor<'static>],
        cache: &mut FieldIdCache,
    ) -> Result<()> {
        for descriptor in info_fields {
            if !descriptor.required {
                continue;
            }

            if descriptor.try_read_field_id(record, cache).is_some() {
                continue;
            }

            let field_name = String::from_utf8_lossy(descriptor.id.as_ref());
            return Err(format!(
                "{} is missing required INFO field {}",
                reader_file_path, field_name
            ));
        }

        Ok(())
    }

    fn read_info_field(
        record: &Record,
        descriptor: &FieldDescriptor<'static>,
        cache: &mut FieldIdCache,
    ) -> Result<Option<InfoFieldValue>> {
        let Some(field_id) = descriptor.try_read_field_id(record, cache) else {
            return Ok(None);
        };

        match descriptor.vcf_type {
            VcfType::String | VcfType::Character => {
                let values = record.info(field_id).string().map_err(|e| e.to_string())?;
                Ok(values.map(|buffer| {
                    let owned = buffer.iter().map(|value| value.to_vec()).collect();
                    InfoFieldValue::String(owned)
                }))
            }
            VcfType::Integer => {
                let values = record.info(field_id).integer().map_err(|e| e.to_string())?;
                Ok(values.map(|buffer| InfoFieldValue::Integer(buffer.to_vec())))
            }
            VcfType::Float => {
                let values = record.info(field_id).float().map_err(|e| e.to_string())?;
                Ok(values.map(|buffer| InfoFieldValue::Float(buffer.to_vec())))
            }
            VcfType::Flag => {
                let flag = record.info(field_id).flag().map_err(|e| e.to_string())?;
                if flag {
                    Ok(Some(InfoFieldValue::Flag(true)))
                } else {
                    Ok(None)
                }
            }
        }
    }

    fn merge_and_write(
        &mut self,
        collected: CollectedVariantData,
        sample_records: &[Option<Record>],
    ) -> Result<()> {
        let template_record = Self::get_template_record(sample_records);

        let dummy = self.writer.dummy_record_mut();
        dummy.set_rid(template_record.rid());
        dummy.set_pos(template_record.pos());
        dummy.set_qual(MISSING_FLOAT);

        // Temporarily move the descriptors out so we can iterate without cloning
        let info_descriptors = mem::take(&mut self.active_info_fields);
        let set_info_result = (|| -> Result<()> {
            for descriptor in &info_descriptors {
                self.set_info_field_from_all_records(sample_records, descriptor)?;
            }
            Ok(())
        })();
        self.active_info_fields = info_descriptors;
        set_info_result?;

        let (merged_genotypes, merged_alleles) =
            merge_exact(collected.genotypes, collected.alleles.as_slice())?;
        let genotypes_flat = Self::flatten_genotypes(&merged_genotypes);

        self.writer
            .dummy_record_mut()
            .set_alleles(&merged_alleles)
            .map_err(|e| e.to_string())?;

        self.write_format_fields(genotypes_flat)?;

        self.writer.write()?;
        self.writer.clear_dummy_record();

        Ok(())
    }

    fn get_template_record(sample_records: &[Option<Record>]) -> &Record {
        let template_index = sample_records
            .iter()
            .position(|r| r.is_some())
            .expect("at least one record must be present");
        sample_records[template_index].as_ref().unwrap()
    }

    fn set_info_field_from_all_records(
        &mut self,
        sample_records: &[Option<Record>],
        descriptor: &FieldDescriptor<'static>,
    ) -> Result<()> {
        let mut selected_value: Option<InfoFieldValue> = None;
        let mut selected_idx: Option<usize> = None;

        for (reader_idx, record) in sample_records
            .iter()
            .enumerate()
            .filter_map(|(idx, record)| record.as_ref().map(|rec| (idx, rec)))
        {
            let value =
                Self::read_info_field(record, descriptor, &mut self.field_id_caches[reader_idx])?;
            let Some(value) = value else { continue };

            if let Some(existing) = &selected_value {
                if !existing.equals(&value) {
                    let field_name = String::from_utf8_lossy(descriptor.id.as_ref());
                    let first_idx = selected_idx.unwrap();
                    let first_label = self
                        .vcf_readers
                        .readers
                        .get(first_idx)
                        .map(|reader| format!("{} (index {})", reader.file_path, first_idx))
                        .unwrap_or_else(|| format!("reader {}", first_idx));
                    let current_label = self
                        .vcf_readers
                        .readers
                        .get(reader_idx)
                        .map(|reader| format!("{} (index {})", reader.file_path, reader_idx))
                        .unwrap_or_else(|| format!("reader {}", reader_idx));
                    return Err(format!(
                        "Conflicting INFO field {} values: {} reports {}, {} reports {}.",
                        field_name,
                        first_label,
                        existing.describe(),
                        current_label,
                        value.describe(),
                    ));
                }
            } else {
                selected_value = Some(value);
                selected_idx = Some(reader_idx);
            }
        }

        if let Some(value) = selected_value {
            value.write_to_record(&mut self.writer, descriptor)?;
        }

        Ok(())
    }

    fn process_genotypes(
        record: &Record,
        genotypes: &mut Vec<Vec<Vec<GenotypeAllele>>>,
    ) -> Result<()> {
        let gt_field = record.genotypes().map_err(|e| e.to_string())?;
        let file_genotypes: Vec<Vec<GenotypeAllele>> = (0..record.sample_count())
            .map(|i| gt_field.get(i as usize).iter().copied().collect())
            .collect();
        genotypes.push(file_genotypes);
        Ok(())
    }

    fn write_format_fields(&mut self, out_gts: Vec<i32>) -> Result<()> {
        // Note: GT is special and written separately from field data registry
        self.writer
            .dummy_record_mut()
            .push_format_integer(b"GT", &out_gts)
            .map_err(|e| e.to_string())?;

        self.field_registry.write_format_fields(&mut self.writer)
    }

    fn flatten_genotypes(out_gts: &[Vec<Vec<GenotypeAllele>>]) -> Vec<i32> {
        let total_capacity: usize = out_gts
            .iter()
            .flat_map(|file_gts| file_gts.iter())
            .map(|sample_gt| {
                if sample_gt.len() == 1 {
                    2
                } else {
                    sample_gt.len()
                }
            })
            .sum();
        let mut gts_new = Vec::with_capacity(total_capacity);
        for file_gts in out_gts {
            for sample_gt in file_gts {
                if sample_gt.len() == 1 {
                    gts_new.push(i32::from(sample_gt[0]));
                    gts_new.push(VECTOR_END_INTEGER);
                } else {
                    gts_new.extend(sample_gt.iter().map(|gt| i32::from(*gt)));
                }
            }
        }
        gts_new
    }

    pub fn write_index(self) -> Result<()> {
        self.writer.write_index()
    }
}

#[derive(Clone, Debug)]
enum InfoFieldValue {
    String(Vec<Vec<u8>>),
    Integer(Vec<i32>),
    Float(Vec<f32>),
    Flag(bool),
}

impl InfoFieldValue {
    fn equals(&self, other: &Self) -> bool {
        match (self, other) {
            (InfoFieldValue::String(a), InfoFieldValue::String(b)) => a == b,
            (InfoFieldValue::Integer(a), InfoFieldValue::Integer(b)) => a == b,
            (InfoFieldValue::Float(a), InfoFieldValue::Float(b)) => {
                a.len() == b.len()
                    && a.iter()
                        .zip(b.iter())
                        .all(|(x, y)| x.to_bits() == y.to_bits())
            }
            (InfoFieldValue::Flag(a), InfoFieldValue::Flag(b)) => a == b,
            _ => false,
        }
    }

    fn describe(&self) -> String {
        match self {
            InfoFieldValue::String(values) => {
                let parts: Vec<String> = values
                    .iter()
                    .map(|value| String::from_utf8_lossy(value).into_owned())
                    .collect();
                format!("[{}]", parts.join(","))
            }
            InfoFieldValue::Integer(values) => format!("{:?}", values),
            InfoFieldValue::Float(values) => format!("{:?}", values),
            InfoFieldValue::Flag(flag) => flag.to_string(),
        }
    }

    fn write_to_record(
        &self,
        writer: &mut VcfWriter,
        descriptor: &FieldDescriptor<'static>,
    ) -> Result<()> {
        let field_id = descriptor.id.as_ref();
        match (descriptor.vcf_type, self) {
            (VcfType::String | VcfType::Character, InfoFieldValue::String(values)) => {
                let borrowed: Vec<&[u8]> = values.iter().map(|value| value.as_slice()).collect();
                writer
                    .dummy_record_mut()
                    .push_info_string(field_id, &borrowed)
                    .map_err(|e| e.to_string())?
            }
            (VcfType::Integer, InfoFieldValue::Integer(values)) => writer
                .dummy_record_mut()
                .push_info_integer(field_id, values)
                .map_err(|e| e.to_string())?,
            (VcfType::Float, InfoFieldValue::Float(values)) => writer
                .dummy_record_mut()
                .push_info_float(field_id, values)
                .map_err(|e| e.to_string())?,
            (VcfType::Flag, InfoFieldValue::Flag(true)) => writer
                .dummy_record_mut()
                .push_info_flag(field_id)
                .map_err(|e| e.to_string())?,
            (VcfType::Flag, InfoFieldValue::Flag(false)) => {}
            _ => {
                return Err(format!(
                    "Type mismatch when writing INFO field {}",
                    String::from_utf8_lossy(field_id)
                ))
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::merge::field_registry::{FieldDataRegistry, FieldNumber};
    use crate::utils::test_util::{
        write_trgt_vcf, write_trgt_vcf_with_records, TestVcfBuilder, TestVcfRecord,
    };
    use rust_htslib::bcf::{index, Read, Reader};
    use std::collections::BinaryHeap;
    use std::path::Path;
    use tempfile::{tempdir, NamedTempFile};

    type RecordSummary = (String, i64, String, Vec<Vec<u8>>);

    fn record_missing_required_format_field() -> Record {
        // Note: ALLR is intentionally missing
        let temp_file = TestVcfBuilder::new()
            .contig("chr1", 1)
            .add_info("TRID", "1", "String", "Test")
            .add_info("END", "1", "Integer", "Test")
            .add_info("MOTIFS", "1", "String", "Test")
            .add_info("STRUC", "1", "String", "Test")
            .add_format("AL", ".", "Integer", "Allele length")
            .add_format("SD", ".", "Integer", "Spanning reads")
            .add_format("AP", ".", "Float", "Allele purity")
            .sample("Sample1")
            .sample("Sample2")
            .record(
                TestVcfRecord::new()
                    .pos(999)
                    .alleles(&["A", "AT", "ATT"])
                    .info_string("TRID", "test1")
                    .info_integer("END", &[1020])
                    .info_string("MOTIFS", "T")
                    .info_string("STRUC", "(T)n")
                    .genotype(&[
                        GenotypeAllele::Unphased(0),
                        GenotypeAllele::Unphased(1),
                        GenotypeAllele::Unphased(1),
                        GenotypeAllele::Unphased(2),
                    ])
                    .format_integer("AL", &[10, 15, 15, 20])
                    .format_integer("SD", &[5, 8, 8, 10])
                    .format_float("AP", &[0.95, 0.98, 0.97, 0.99]),
            )
            .build();

        let mut reader = Reader::from_path(temp_file.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();
        record
    }

    fn record_missing_required_info_field() -> Record {
        let temp_file = TestVcfBuilder::new()
            .contig("chr1", 1)
            .with_trgt_defaults()
            .with_trgt_format_defaults()
            .sample("Sample1")
            .record(
                TestVcfRecord::new()
                    .pos(999)
                    .alleles(&["A", "AT", "ATT"])
                    .info_string("TRID", "test1")
                    .info_integer("END", &[1020])
                    // Note: MOTIFS intentionally omitted here to test error handling
                    .info_string("STRUC", "(T)n")
                    .genotype(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
                    .format_integer("AL", &[10, 15])
                    .format_string("ALLR", &["8-12", "13-17"])
                    .format_integer("SD", &[5, 8])
                    .format_string("MC", &["10", "15"])
                    .format_string("MS", &["0-10", "0-15"])
                    .format_float("AP", &[0.95, 0.98])
                    .format_float("AM", &[0.5, 0.6]),
            )
            .build();

        let mut reader = Reader::from_path(temp_file.path()).unwrap();
        let mut record = reader.empty_record();
        let _ = reader.read(&mut record).unwrap();
        record
    }

    fn create_mock_processor() -> VcfProcessor {
        let header = TestVcfBuilder::new()
            .contig("chr1", 1)
            .with_trgt_defaults()
            .with_trgt_format_defaults()
            .sample("Sample1")
            .sample("Sample2")
            .build_header();

        let vcf_readers = VcfReaders::empty(VcfReadMode::Indexed);
        let temp_output = NamedTempFile::new().unwrap();
        let writer =
            VcfWriter::new(&header, None, Some(&temp_output.path().to_path_buf()), None).unwrap();

        let active_format_fields: Vec<FieldDescriptor<'static>> = FORMAT_FIELDS
            .iter()
            .filter(|d| !d.special_handling)
            .cloned()
            .collect();
        let field_registry = FieldDataRegistry::from_descriptors(&active_format_fields);

        // Tests may pass sample_records with more entries than actual readers so pre-allocate enough caches for typical test cases
        const TEST_CACHE_COUNT: usize = 4;
        VcfProcessor {
            writer,
            vcf_readers,
            genome_reader: None,
            contig_order: vec!["chr1".to_string()],
            contig_to_idx: HashMap::from([("chr1".to_string(), 0)]),
            read_mode: VcfReadMode::Indexed,
            skip_n: 0,
            process_n: usize::MAX,
            needs_padding: false,
            quit_on_error: false,
            active_format_fields,
            active_info_fields: INFO_FIELDS.to_vec(),
            field_registry,
            field_id_caches: (0..TEST_CACHE_COUNT).map(|_| FieldIdCache::new()).collect(),
        }
    }

    fn missing_required_format_field_error(processor: &VcfProcessor) -> String {
        let record = record_missing_required_format_field();
        let mut field_data = FieldDataRegistry::from_descriptors(&processor.active_format_fields);
        let mut cache = FieldIdCache::new();

        field_data
            .append_record(&record, &processor.active_format_fields, &mut cache)
            .expect_err("missing required ALLR field should error")
    }

    fn missing_required_info_field_error(processor: &mut VcfProcessor) -> String {
        let sample_records = vec![Some(record_missing_required_info_field())];
        match processor.collect_variant_data(&sample_records) {
            Ok(_) => panic!("collect_variant_data should fail when MOTIFS is missing"),
            Err(err) => err,
        }
    }

    fn merge_args_for_test(vcfs: Vec<PathBuf>, output: PathBuf, no_index: bool) -> MergeArgs {
        MergeArgs {
            vcfs: Some(vcfs),
            vcf_list: None,
            genome_src: None,
            output: Some(output),
            output_type: None,
            skip_n: None,
            process_n: None,
            print_header: false,
            force_single: false,
            no_version: true,
            quit_on_error: true,
            contigs: None,
            threads: 1,
            no_index: no_index.into(),
            write_index: false,
        }
    }

    fn summarize_vcf(path: &Path) -> Result<(usize, Vec<RecordSummary>)> {
        let mut reader = Reader::from_path(path).map_err(|e| e.to_string())?;
        let header = reader.header().clone();
        let sample_count = header.sample_count() as usize;
        let mut rows = Vec::new();

        for record in reader.records() {
            let record = record.map_err(|e| e.to_string())?;
            let rid = record
                .rid()
                .ok_or_else(|| "Record missing RID".to_string())?;
            let contig = header.rid2name(rid).map_err(|e| e.to_string())?;
            let trid = record
                .info(b"TRID")
                .string()
                .map_err(|e| e.to_string())?
                .ok_or_else(|| "Missing TRID".to_string())?;
            let trid_str = String::from_utf8_lossy(trid[0]).into_owned();
            let alleles = record.alleles().iter().map(|a| a.to_vec()).collect();

            rows.push((
                String::from_utf8_lossy(contig).into_owned(),
                record.pos(),
                trid_str,
                alleles,
            ));
        }

        Ok((sample_count, rows))
    }

    fn make_test_records(trids: &[&str]) -> Vec<Record> {
        let mut builder = TestVcfBuilder::new()
            .with_trgt_defaults()
            .contig("chr1", 10000);

        for trid in trids {
            builder = builder.record(
                TestVcfRecord::with_trgt_info(trid.as_bytes(), 1020, b"T", b"(T)n")
                    .pos(999)
                    .alleles(&["A", "AT"]),
            );
        }

        let temp_file = builder.build();
        let mut reader = Reader::from_path(temp_file.path()).unwrap();
        let mut records = Vec::new();
        loop {
            let mut record = reader.empty_record();
            match reader.read(&mut record) {
                Some(Ok(())) => records.push(record),
                _ => break,
            }
        }
        records
    }

    #[test]
    fn test_vcf_record_wrapper_heap() {
        let temp_file = TestVcfBuilder::new().build();
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

    #[test]
    fn test_merge_matches_indexed_for_simple_inputs() -> Result<()> {
        let dir = tempdir().unwrap();
        let vcf_a = dir.path().join("a.vcf.gz");
        let vcf_b = dir.path().join("b.vcf.gz");

        write_trgt_vcf(&vcf_a, "S1", "repeat1", 10)?;
        write_trgt_vcf(&vcf_b, "S2", "repeat1", 10)?;

        index::build(&vcf_a, None::<&PathBuf>, 1, index::Type::Tbx).map_err(|e| e.msg)?;
        index::build(&vcf_b, None::<&PathBuf>, 1, index::Type::Tbx).map_err(|e| e.msg)?;

        let output_indexed = dir.path().join("indexed.vcf");
        let output_stream = dir.path().join("stream.vcf");

        let vcfs = vec![vcf_a.clone(), vcf_b.clone()];
        let args_indexed = merge_args_for_test(vcfs.clone(), output_indexed.clone(), false);
        let args_stream = merge_args_for_test(vcfs.clone(), output_stream.clone(), true);

        let mut processor_indexed = VcfProcessor::new(&args_indexed, vcfs.clone())?;
        processor_indexed.merge_variants()?;

        let mut processor_stream = VcfProcessor::new(&args_stream, vcfs)?;
        processor_stream.merge_variants()?;

        drop(processor_indexed);
        drop(processor_stream);

        let (indexed_samples, indexed_records) = summarize_vcf(&output_indexed)?;
        let (stream_samples, stream_records) = summarize_vcf(&output_stream)?;

        assert_eq!(indexed_samples, 2);
        assert_eq!(indexed_samples, stream_samples);
        assert_eq!(indexed_records, stream_records);

        Ok(())
    }

    #[test]
    fn test_duplicate_records_within_reader() -> Result<()> {
        let dir = tempdir().unwrap();
        let vcf_path = dir.path().join("dupe.vcf");

        write_trgt_vcf_with_records(
            vec![
                TestVcfRecord::default_trgt("repeatX", 1000, 1050),
                TestVcfRecord::default_trgt("repeatX", 1000, 1050),
            ],
            &vcf_path,
            "S1",
        )?;

        let output = dir.path().join("out.vcf");
        let mut args = merge_args_for_test(vec![vcf_path.clone()], output, true);
        args.force_single = true;

        let mut processor = VcfProcessor::new(&args, vec![vcf_path])?;
        let err = processor
            .merge_variants()
            .expect_err("duplicate records at same site should error in streaming mode");

        assert!(
            err.contains("Multiple records detected"),
            "error should mention duplicate detection, got: {}",
            err
        );

        Ok(())
    }

    #[test]
    fn test_streaming_merges_all_records_at_same_position() -> Result<()> {
        let dir = tempdir().unwrap();
        let vcf_a = dir.path().join("a_multi.vcf.gz");
        let vcf_b = dir.path().join("b_multi.vcf.gz");

        write_trgt_vcf_with_records(
            vec![
                TestVcfRecord::default_trgt("repeatA", 1000, 1050),
                TestVcfRecord::default_trgt("repeatB", 1000, 1050),
            ],
            &vcf_a,
            "S1",
        )?;

        write_trgt_vcf_with_records(
            vec![TestVcfRecord::default_trgt("repeatB", 1000, 1050)],
            &vcf_b,
            "S2",
        )?;

        index::build(&vcf_a, None::<&PathBuf>, 1, index::Type::Tbx).map_err(|e| e.msg)?;
        index::build(&vcf_b, None::<&PathBuf>, 1, index::Type::Tbx).map_err(|e| e.msg)?;

        let output_indexed = dir.path().join("indexed_multi.vcf");
        let output_stream = dir.path().join("stream_multi.vcf");

        let vcfs = vec![vcf_a.clone(), vcf_b.clone()];
        let args_indexed = merge_args_for_test(vcfs.clone(), output_indexed.clone(), false);
        let args_stream = merge_args_for_test(vcfs.clone(), output_stream.clone(), true);

        let mut processor_indexed = VcfProcessor::new(&args_indexed, vcfs.clone())?;
        processor_indexed.merge_variants()?;

        let mut processor_stream = VcfProcessor::new(&args_stream, vcfs)?;
        processor_stream.merge_variants()?;

        drop(processor_indexed);
        drop(processor_stream);

        let (indexed_samples, indexed_records) = summarize_vcf(&output_indexed)?;
        let (stream_samples, stream_records) = summarize_vcf(&output_stream)?;

        assert_eq!(indexed_samples, 2);
        assert_eq!(indexed_samples, stream_samples);
        assert_eq!(indexed_records, stream_records);

        Ok(())
    }

    #[test]
    fn test_streaming_merge_errors_on_decreasing_positions() -> Result<()> {
        let dir = tempdir().unwrap();
        let vcf_path = dir.path().join("unsorted.vcf");

        write_trgt_vcf_with_records(
            vec![
                TestVcfRecord::default_trgt("repeat1", 1000, 1050),
                TestVcfRecord::default_trgt("repeat1", 500, 550),
            ],
            &vcf_path,
            "S1",
        )?;

        let output = dir.path().join("out.vcf");
        let mut args = merge_args_for_test(vec![vcf_path.clone()], output, true);
        args.force_single = true;

        let mut processor = VcfProcessor::new(&args, vec![vcf_path])?;
        let err = processor
            .merge_variants()
            .expect_err("streaming merge should fail on decreasing positions");

        assert!(
            err.contains("not sorted correctly within contig"),
            "error should mention contig sort failure, got: {}",
            err
        );

        Ok(())
    }

    #[test]
    fn test_missing_required_fields_error() {
        let mut processor = create_mock_processor();
        let err = missing_required_format_field_error(&processor);
        assert!(
            err.contains("ALLR"),
            "error should mention missing ALLR, got: {}",
            err
        );

        let err = missing_required_info_field_error(&mut processor);
        assert!(
            err.contains("MOTIFS"),
            "Error message should mention MOTIFS, got: {}",
            err
        );
    }

    #[test]
    fn test_merge_and_write_restores_info_fields_after_error() {
        let mut processor = create_mock_processor();
        assert!(
            !processor.active_info_fields.is_empty(),
            "precondition: info descriptors should be present"
        );

        let mut records = make_test_records(&["repeat1", "repeat2"]);
        let sample_records = vec![Some(records.remove(0)), Some(records.remove(0))];
        let collected = CollectedVariantData {
            genotypes: vec![
                vec![vec![GenotypeAllele::Unphased(0)]],
                vec![vec![GenotypeAllele::Unphased(0)]],
            ],
            alleles: vec![
                vec![b"A".as_ref(), b"AT".as_ref()],
                vec![b"A".as_ref(), b"AT".as_ref()],
            ],
        };

        let err = processor
            .merge_and_write(collected, sample_records.as_slice())
            .expect_err("conflicting TRIDs should surface as an error");

        assert!(
            err.contains("Conflicting INFO field TRID"),
            "error message should mention TRID conflict, got: {}",
            err
        );
        assert_eq!(
            processor.active_info_fields.len(),
            INFO_FIELDS.len(),
            "INFO descriptors should be restored even after an error"
        );
        assert!(
            processor
                .writer
                .dummy_record()
                .info(b"TRID")
                .string()
                .expect("info lookup should succeed")
                .is_none(),
            "dummy record should not retain TRID after a failed merge"
        );
    }

    #[test]
    fn test_set_info_field_prefers_first_non_missing_value() {
        // Record without TRID (first), followed by two records with TRID
        let temp_file = TestVcfBuilder::new()
            .contig("chr1", 10000)
            .with_trgt_defaults()
            .record(
                TestVcfRecord::new()
                    .pos(999)
                    .alleles(&["A", "AT"])
                    .info_integer("END", &[1020])
                    .info_string("MOTIFS", "T")
                    .info_string("STRUC", "(T)n"),
            )
            .record(
                TestVcfRecord::new()
                    .pos(999)
                    .alleles(&["A", "AT"])
                    .info_string("TRID", "repeatX")
                    .info_integer("END", &[1020])
                    .info_string("MOTIFS", "T")
                    .info_string("STRUC", "(T)n"),
            )
            .record(
                TestVcfRecord::new()
                    .pos(999)
                    .alleles(&["A", "AT"])
                    .info_string("TRID", "repeatX")
                    .info_integer("END", &[1020])
                    .info_string("MOTIFS", "T")
                    .info_string("STRUC", "(T)n"),
            )
            .build();

        let mut reader = Reader::from_path(temp_file.path()).unwrap();
        let sample_records: Vec<Option<Record>> = (0..3)
            .map(|_| {
                let mut record = reader.empty_record();
                let _ = reader.read(&mut record);
                Some(record)
            })
            .collect();

        let descriptor = INFO_FIELDS
            .iter()
            .find(|d| d.id.as_ref() == b"TRID")
            .unwrap();
        let mut processor = create_mock_processor();
        processor
            .set_info_field_from_all_records(sample_records.as_slice(), descriptor)
            .expect("should select first available TRID");

        let values = processor
            .writer
            .dummy_record()
            .info(b"TRID")
            .string()
            .expect("info lookup should succeed")
            .expect("TRID value should be present");
        assert_eq!(values.len(), 1);
        assert_eq!(values[0], b"repeatX");
    }

    #[test]
    fn test_custom_info_field_optional_across_vcfs() {
        let temp_with_custom = TestVcfBuilder::new()
            .contig("chr1", 10000)
            .add_info("CUSTOM", "1", "String", "Custom")
            .record(
                TestVcfRecord::new()
                    .pos(999)
                    .alleles(&["A", "AT"])
                    .info_string("CUSTOM", "value"),
            )
            .build();

        let temp_without_custom = TestVcfBuilder::new()
            .contig("chr1", 10000)
            .record(TestVcfRecord::new().pos(999).alleles(&["A", "AT"]))
            .build();

        let mut reader_with_custom = Reader::from_path(temp_with_custom.path()).unwrap();
        let mut record_with_custom = reader_with_custom.empty_record();
        let _ = reader_with_custom.read(&mut record_with_custom);

        let mut reader_without_custom = Reader::from_path(temp_without_custom.path()).unwrap();
        let mut record_without_custom = reader_without_custom.empty_record();
        let _ = reader_without_custom.read(&mut record_without_custom);

        let sample_records = vec![Some(record_with_custom), Some(record_without_custom)];

        let custom_descriptor = FieldDescriptor::new(
            b"CUSTOM",
            FieldCategory::Info,
            VcfType::String,
            FieldNumber::Fixed(1),
            "Custom optional field",
        );

        let header = TestVcfBuilder::new()
            .contig("chr1", 10000)
            .add_info("CUSTOM", "1", "String", "Custom optional field")
            .sample("S1")
            .build_header();

        let temp_output = NamedTempFile::new().unwrap();
        let writer =
            VcfWriter::new(&header, None, Some(&temp_output.path().to_path_buf()), None).unwrap();
        let active_format_fields: Vec<FieldDescriptor<'static>> = Vec::new();
        let field_registry = FieldDataRegistry::from_descriptors(&active_format_fields);

        let mut processor = VcfProcessor {
            writer,
            vcf_readers: VcfReaders::empty(VcfReadMode::Indexed),
            genome_reader: None,
            contig_order: vec!["chr1".to_string()],
            contig_to_idx: HashMap::from([("chr1".to_string(), 0)]),
            read_mode: VcfReadMode::Indexed,
            skip_n: 0,
            process_n: usize::MAX,
            needs_padding: false,
            quit_on_error: false,
            active_format_fields,
            active_info_fields: vec![custom_descriptor.clone()],
            field_registry,
            field_id_caches: vec![FieldIdCache::new(), FieldIdCache::new()],
        };

        processor
            .set_info_field_from_all_records(&sample_records, &custom_descriptor)
            .expect("custom INFO field should be optional");

        let values = processor
            .writer
            .dummy_record()
            .info(b"CUSTOM")
            .string()
            .expect("custom field should be present")
            .expect("custom field value should exist");
        assert_eq!(values.len(), 1);
        assert_eq!(values[0], b"value");
    }

    #[test]
    fn test_group_records_by_trid_handles_multiple_trids() {
        let processor = create_mock_processor();
        let mut per_reader_records = vec![
            make_test_records(&["repeatB", "repeatA"]),
            make_test_records(&["repeatA"]),
        ];

        let grouped = processor
            .group_records_by_trid(&mut per_reader_records, "chr1", 1000)
            .expect("grouping should succeed");

        let trids: Vec<String> = grouped
            .keys()
            .map(|k| String::from_utf8_lossy(k).into_owned())
            .collect();
        assert_eq!(trids, vec!["repeatA".to_string(), "repeatB".to_string()]);

        assert_eq!(
            grouped.get("repeatA".as_bytes()).unwrap().len(),
            2,
            "repeatA should gather records from both readers"
        );
        assert_eq!(
            grouped.get("repeatB".as_bytes()).unwrap().len(),
            1,
            "repeatB should only exist in reader 0"
        );
    }

    #[test]
    fn test_load_group_into_sample_records_detects_duplicates() {
        let processor = create_mock_processor();
        let mut reader_records = make_test_records(&["repeatA", "repeatA"]);

        let mut group_records = Vec::new();
        for record in reader_records.drain(..) {
            group_records.push((0usize, record));
        }

        let mut sample_records: Vec<Option<Record>> = vec![None];
        let err = processor
            .load_group_into_sample_records(
                &mut sample_records,
                group_records,
                "chr1",
                1000,
                b"repeatA",
            )
            .expect_err("duplicate reader entries should error");
        assert!(
            err.contains("repeatA"),
            "error should reference TRID, got {err}"
        );
    }
}
