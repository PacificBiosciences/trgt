use super::tpool::HtsThreadPool;
use crate::utils::{OutputType, Result};
use rust_htslib::bcf;
use std::{path::PathBuf, sync::Arc};

pub struct VcfWriter {
    // writer must drop BEFORE _thread_pool to ensure pool is alive during flush
    writer: bcf::Writer,
    _thread_pool: Option<Arc<HtsThreadPool>>,
    dummy_record: bcf::Record,
    output_type: OutputType,
    output_path: Option<PathBuf>,
}

impl VcfWriter {
    pub fn new(
        header: &bcf::Header,
        output_type: Option<&OutputType>,
        output: Option<&PathBuf>,
        thread_pool: Option<Arc<HtsThreadPool>>,
    ) -> Result<Self> {
        let output_type = output_type
            .cloned()
            .or_else(|| output.map(|path| OutputType::infer_from_path(path)))
            .unwrap_or(OutputType::Vcf {
                is_uncompressed: true,
                level: None,
            });
        log::debug!("Creating {} writer", output_type);

        let (is_uncompressed, format) = output_type.writer_options();
        let mut writer = match output {
            Some(path) => bcf::Writer::from_path(path, header, is_uncompressed, format),
            None => bcf::Writer::from_stdout(header, is_uncompressed, format),
        }
        .map_err(|e| format!("Failed to create writer: {}", e))?;

        if !is_uncompressed {
            if let Some(ref pool) = thread_pool {
                // SAFETY: Writer is created and owned here and the pool lives alongside it
                // and attach_to_writer only forwards the inner htsFile pointer
                unsafe {
                    pool.attach_to_writer(&mut writer)?;
                }
            }
        }

        let dummy_record = writer.empty_record();
        let output_path = output.cloned();

        Ok(VcfWriter {
            writer,
            _thread_pool: thread_pool,
            dummy_record,
            output_type,
            output_path,
        })
    }

    pub fn dummy_record(&self) -> &bcf::Record {
        &self.dummy_record
    }

    pub fn dummy_record_mut(&mut self) -> &mut bcf::Record {
        &mut self.dummy_record
    }

    pub fn write(&mut self) -> Result<()> {
        self.writer
            .write(&self.dummy_record)
            .map_err(|e| e.to_string())
    }

    pub fn clear_dummy_record(&mut self) {
        self.dummy_record.clear();
    }

    fn into_parts(self) -> (bcf::Writer, OutputType, Option<PathBuf>) {
        let VcfWriter {
            writer,
            output_type,
            output_path,
            ..
        } = self;
        (writer, output_type, output_path)
    }

    pub fn write_index(self) -> Result<()> {
        let (writer, output_type, output_path) = self.into_parts();

        let Some(output_path) = output_path else {
            log::debug!("Skipping index creation because --output was not provided");
            return Ok(());
        };

        // Drop writer so the BGZF footer is flushed before indexing the file.
        drop(writer);

        let index_type = match output_type {
            OutputType::Vcf {
                is_uncompressed, ..
            } if !is_uncompressed => Some((bcf::index::Type::Tbx, "tabix")),
            OutputType::Bcf {
                is_uncompressed, ..
            } if !is_uncompressed => Some((bcf::index::Type::Csi(14), "CSI")),
            _ => None,
        };

        let Some((index_type, index_label)) = index_type else {
            log::debug!("Skipping index creation for {} output types", output_type);
            return Ok(());
        };

        log::debug!(
            "Building {} index for {}",
            index_label,
            output_path.display()
        );

        bcf::index::build(&output_path, None, 1, index_type)
            .map_err(|e| format!("Failed to build index for {}: {}", output_path.display(), e))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::test_util::TestVcfBuilder;
    use rust_htslib::bcf::record::GenotypeAllele;
    use std::path::Path;
    use tempfile::{tempdir, NamedTempFile};

    fn populate_dummy_record(writer: &mut VcfWriter) {
        let dummy = writer.dummy_record_mut();
        dummy.set_rid(Some(0));
        dummy.set_pos(100);
        dummy.set_qual(42.0);
        dummy
            .set_alleles(&[b"A", b"T"])
            .expect("alleles should set up");
        dummy
            .push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
            .expect("genotype push succeeds");
    }

    fn expected_index_path(output_path: &Path, suffix: &str) -> PathBuf {
        output_path.parent().unwrap().join(format!(
            "{}.{suffix}",
            output_path.file_name().unwrap().to_string_lossy()
        ))
    }

    #[test]
    fn dummy_record_helpers_write_and_clear() {
        let header = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .sample("S1")
            .build_header();

        let temp_output = NamedTempFile::new().unwrap();
        let mut writer =
            VcfWriter::new(&header, None, Some(&temp_output.path().to_path_buf()), None).unwrap();

        populate_dummy_record(&mut writer);

        writer.write().expect("write succeeds");
        writer.clear_dummy_record();
        assert_eq!(writer.dummy_record().alleles().len(), 0);
    }

    #[test]
    fn write_index_builds_tabix_for_compressed_vcf() {
        let header = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .sample("S1")
            .build_header();
        let temp_dir = tempdir().unwrap();
        let output_path = temp_dir.path().join("out.vcf.gz");

        let mut writer = VcfWriter::new(&header, None, Some(&output_path), None).unwrap();
        populate_dummy_record(&mut writer);
        writer.write().expect("vcf write succeeds");

        let expected_index = expected_index_path(&output_path, "tbi");
        assert!(
            !expected_index.exists(),
            "index should not exist before write_index"
        );

        writer.write_index().expect("index build should succeed");
        assert!(
            expected_index.exists(),
            "tabix index should exist after write_index"
        );
    }

    #[test]
    fn write_index_builds_csi_for_compressed_bcf() {
        let header = TestVcfBuilder::new()
            .contig("chr1", 1000)
            .sample("S1")
            .build_header();
        let temp_dir = tempdir().unwrap();
        let output_path = temp_dir.path().join("out.bcf.gz");
        let output_type = OutputType::Bcf {
            is_uncompressed: false,
            level: None,
        };

        let mut writer =
            VcfWriter::new(&header, Some(&output_type), Some(&output_path), None).unwrap();
        populate_dummy_record(&mut writer);
        writer.write().expect("bcf write succeeds");

        let expected_index = expected_index_path(&output_path, "csi");
        assert!(
            !expected_index.exists(),
            "index should not exist before write_index"
        );

        writer.write_index().expect("index build should succeed");
        assert!(
            expected_index.exists(),
            "csi index should exist after write_index"
        );
    }
}
