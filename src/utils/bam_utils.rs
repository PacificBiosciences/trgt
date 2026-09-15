use crate::utils::{open_bam_reader, BamSources, InputSource, Result};
use rust_htslib::bam::{self, Read};
use std::collections::{BTreeSet, HashSet};

pub fn get_bam_header(reads_src: &InputSource) -> Result<bam::Header> {
    let bam = open_bam_reader(reads_src, 1)?;
    Ok(bam::Header::from_template(bam.header()))
}

pub fn has_reference_targets(bam_header: &bam::Header) -> bool {
    bam::HeaderView::from_header(bam_header).target_count() > 0
}

pub fn get_sample_name(reads_src: &InputSource, bam_header: &bam::Header) -> Result<String> {
    let header_hashmap = bam_header.to_hashmap();
    let mut sample_names = HashSet::new();

    if let Some(rg_fields) = header_hashmap.get("RG") {
        for rg_field in rg_fields {
            if let Some(sample_name) = rg_field.get("SM") {
                sample_names.insert(sample_name.to_owned());
            }
        }
    }

    match sample_names.len() {
        1 => {
            return Ok(sample_names
                .into_iter()
                .next()
                .expect("sample_names has length 1"));
        }
        0 => log::warn!("No sample names found. Falling back to BAM path stem to infer sample name."),
        _ => log::warn!(
            "Multiple sample names found (n = {sample_count}). Falling back to BAM path stem to infer sample name.",
            sample_count = sample_names.len()
        ),
    };

    let sample = reads_src
        .file_stem()
        .ok_or_else(|| "Invalid BAM source; cannot infer sample name from path/URL".to_string())?;

    Ok(sample)
}

pub fn validate_bam_sources(
    bam_sources: &BamSources,
    sample_name_override: Option<&str>,
) -> Result<(bam::Header, String)> {
    reject_duplicate_bam_sources(bam_sources.as_slice())?;

    let mut headers = Vec::with_capacity(bam_sources.as_slice().len());
    for source in bam_sources.as_slice() {
        let reader = open_bam_reader(source, 1)
            .map_err(|e| format!("Failed to validate BAM source {source}: {e}"))?;
        let header = bam::Header::from_template(reader.header());
        if !has_reference_targets(&header) {
            return Err(format!(
                "BAM source {source} is not mapped: header contains no @SQ reference records"
            ));
        }

        if let Some(primary_header) = headers.first() {
            validate_reference_dictionary(bam_sources.primary(), source, primary_header, &header)?;
        }

        headers.push(header);
    }

    let sample_name = resolve_sample_name(bam_sources.as_slice(), &headers, sample_name_override)?;
    let primary_header = headers
        .into_iter()
        .next()
        .expect("BAM sources always contain a primary source");

    Ok((primary_header, sample_name))
}

fn reject_duplicate_bam_sources(bam_sources: &[InputSource]) -> Result<()> {
    let mut seen = HashSet::with_capacity(bam_sources.len());
    for source in bam_sources {
        if !seen.insert(bam_source_identity(source)) {
            return Err(format!("Duplicate BAM source: {source}"));
        }
    }
    Ok(())
}

fn bam_source_identity(source: &InputSource) -> String {
    match source {
        InputSource::Local(path) => {
            let path = std::fs::canonicalize(path).unwrap_or_else(|_| path.clone());
            format!("local:{}", path.display())
        }
        InputSource::Remote(remote) => format!("{}:{}", remote.scheme_str(), remote.url()),
    }
}

fn validate_reference_dictionary(
    primary_source: &InputSource,
    source: &InputSource,
    primary_header: &bam::Header,
    header: &bam::Header,
) -> Result<()> {
    if reference_dictionary(primary_header) != reference_dictionary(header) {
        return Err(format!(
            "BAM source {source} has a reference sequence dictionary that differs from primary BAM source {primary_source}; contig names, order, and lengths must match"
        ));
    }
    Ok(())
}

fn reference_dictionary(header: &bam::Header) -> Vec<(Vec<u8>, u64)> {
    let header_view = bam::HeaderView::from_header(header);
    (0..header_view.target_count())
        .map(|tid| {
            (
                header_view.tid2name(tid).to_vec(),
                header_view
                    .target_len(tid)
                    .expect("target ids below target_count have lengths"),
            )
        })
        .collect()
}

fn resolve_sample_name(
    bam_sources: &[InputSource],
    headers: &[bam::Header],
    sample_name_override: Option<&str>,
) -> Result<String> {
    if bam_sources.len() == 1 {
        return match sample_name_override {
            Some(sample_name) => Ok(sample_name.to_owned()),
            None => get_sample_name(&bam_sources[0], &headers[0]),
        };
    }

    let sample_names = headers
        .iter()
        .flat_map(sample_names_from_header)
        .collect::<BTreeSet<_>>();

    match sample_names.len() {
        0 => sample_name_override.map(str::to_owned).ok_or_else(|| {
            "No sample names found across BAM sources; provide --sample-name for multi-source genotype runs".to_string()
        }),
        1 => Ok(sample_name_override
            .map(str::to_owned)
            .unwrap_or_else(|| sample_names.into_iter().next().expect("sample_names has length 1"))),
        _ => {
            let sample_names = sample_names.into_iter().collect::<Vec<_>>().join(", ");
            Err(format!(
                "Conflicting sample names found across BAM sources: {sample_names}"
            ))
        }
    }
}

fn sample_names_from_header(header: &bam::Header) -> BTreeSet<String> {
    header
        .to_hashmap()
        .get("RG")
        .into_iter()
        .flat_map(|rg_fields| rg_fields.iter())
        .filter_map(|rg_field| rg_field.get("SM"))
        .filter(|sample_name| !sample_name.is_empty())
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::{
        test_util::{
            create_indexed_test_bam as indexed_bam, create_unindexed_test_bam as unindexed_bam,
            TestBam,
        },
        BamSources,
    };
    use rust_htslib::bam;
    use std::path::PathBuf;

    fn bam_sources(primary: &TestBam, supplemental: Option<&TestBam>) -> BamSources {
        BamSources::new(primary.source(), supplemental.map(TestBam::source))
    }

    fn validate(sources: &BamSources, sample_name: Option<&str>) -> Result<(bam::Header, String)> {
        validate_bam_sources(sources, sample_name)
    }

    #[test]
    fn duplicate_bam_sources_are_rejected_before_opening() {
        let source = InputSource::Local(PathBuf::from("/missing/duplicate.bam"));
        let sources = BamSources::new(source.clone(), Some(source.clone()));

        let error = validate(&sources, Some("sample")).unwrap_err();

        assert!(error.contains("Duplicate BAM source"));
        assert!(error.contains("/missing/duplicate.bam"));
    }

    #[test]
    fn unreadable_or_unindexed_supplemental_source_reports_source() {
        let primary = indexed_bam("primary.bam", &[("chr1", 1000)], &["sample"]);
        let supplemental = unindexed_bam("failed.bam", &[("chr1", 1000)], &["sample"]);
        let sources = bam_sources(&primary, Some(&supplemental));

        let error = validate(&sources, None).unwrap_err();

        assert!(error.contains("failed.bam"));
        assert!(error.contains("Failed to create BAM reader"));
    }

    #[test]
    fn headers_without_reference_records_are_rejected() {
        let primary = indexed_bam("primary.bam", &[], &["sample"]);
        let sources = bam_sources(&primary, None);

        let error = validate(&sources, None).unwrap_err();

        assert!(error.contains("primary.bam"));
        assert!(error.contains("no @SQ reference records"));
    }

    #[test]
    fn supplemental_reference_dictionary_must_match_primary() {
        for (name, contigs) in [
            ("name", vec![("chr2", 1000)]),
            ("order", vec![("chr2", 2000), ("chr1", 1000)]),
            ("length", vec![("chr1", 999), ("chr2", 2000)]),
        ] {
            let primary = indexed_bam(
                &format!("primary-{name}.bam"),
                &[("chr1", 1000), ("chr2", 2000)],
                &["sample"],
            );
            let supplemental = indexed_bam(&format!("failed-{name}.bam"), &contigs, &["sample"]);
            let sources = bam_sources(&primary, Some(&supplemental));

            let error = validate(&sources, None).unwrap_err();

            assert!(error.contains(&format!("failed-{name}.bam")));
            assert!(error.contains("reference sequence dictionary"));
        }
    }

    #[test]
    fn multi_source_sample_names_resolve_when_metadata_agrees() {
        let primary = indexed_bam("primary.bam", &[("chr1", 1000)], &["sample"]);
        let supplemental = indexed_bam("failed.bam", &[("chr1", 1000)], &["sample"]);
        let sources = bam_sources(&primary, Some(&supplemental));

        let (_, sample_name) = validate(&sources, None).unwrap();

        assert_eq!(sample_name, "sample");
    }

    #[test]
    fn validation_returns_primary_header_as_canonical_header() {
        let primary = indexed_bam("primary.bam", &[("chr1", 1000)], &["sample"]);
        let supplemental = indexed_bam("failed.bam", &[("chr1", 1000)], &["sample", "sample"]);
        let sources = bam_sources(&primary, Some(&supplemental));

        let (primary_header, _) = validate(&sources, None).unwrap();

        assert_eq!(
            primary_header.to_hashmap(),
            get_bam_header(&primary.source()).unwrap().to_hashmap()
        );
    }

    #[test]
    fn multi_source_sample_name_override_is_output_label_only() {
        let primary = indexed_bam("primary.bam", &[("chr1", 1000)], &["sample"]);
        let supplemental = indexed_bam("failed.bam", &[("chr1", 1000)], &["sample"]);
        let sources = bam_sources(&primary, Some(&supplemental));

        let (_, sample_name) = validate(&sources, Some("output-label")).unwrap();

        assert_eq!(sample_name, "output-label");
    }

    #[test]
    fn multi_source_conflicting_sample_metadata_is_rejected_even_with_override() {
        let primary = indexed_bam("primary.bam", &[("chr1", 1000)], &["sample-a"]);
        let supplemental = indexed_bam("failed.bam", &[("chr1", 1000)], &["sample-b"]);
        let sources = bam_sources(&primary, Some(&supplemental));

        let error = validate(&sources, Some("output-label")).unwrap_err();

        assert!(error.contains("Conflicting sample names"));
        assert!(error.contains("sample-a"));
        assert!(error.contains("sample-b"));
    }

    #[test]
    fn multi_source_allows_missing_sample_metadata_when_another_source_has_sample_name() {
        let primary = indexed_bam("primary.bam", &[("chr1", 1000)], &["sample"]);
        let supplemental = indexed_bam("failed.bam", &[("chr1", 1000)], &[]);
        let sources = bam_sources(&primary, Some(&supplemental));

        let (_, sample_name) = validate(&sources, None).unwrap();

        assert_eq!(sample_name, "sample");
    }

    #[test]
    fn multi_source_requires_explicit_sample_name_when_all_metadata_is_missing() {
        let primary = indexed_bam("primary.bam", &[("chr1", 1000)], &[]);
        let supplemental = indexed_bam("failed.bam", &[("chr1", 1000)], &[]);
        let sources = bam_sources(&primary, Some(&supplemental));

        let error = validate(&sources, None).unwrap_err();

        assert!(error.contains("No sample names found"));

        let (_, sample_name) = validate(&sources, Some("explicit-sample")).unwrap();
        assert_eq!(sample_name, "explicit-sample");
    }

    #[test]
    fn single_source_sample_name_behavior_remains_backward_compatible() {
        let without_metadata = indexed_bam("primary.bam", &[("chr1", 1000)], &[]);
        let sources = bam_sources(&without_metadata, None);
        let (_, sample_name) = validate(&sources, None).unwrap();
        assert_eq!(sample_name, without_metadata.source().file_stem().unwrap());

        let with_multiple_samples =
            indexed_bam("multi.bam", &[("chr1", 1000)], &["sample-a", "sample-b"]);
        let sources = bam_sources(&with_multiple_samples, None);
        let (_, sample_name) = validate(&sources, None).unwrap();
        assert_eq!(
            sample_name,
            with_multiple_samples.source().file_stem().unwrap()
        );

        let (_, sample_name) = validate(&sources, Some("output-label")).unwrap();
        assert_eq!(sample_name, "output-label");
    }
}
