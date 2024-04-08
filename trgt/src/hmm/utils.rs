use super::Span;

pub fn count_motifs(motifs: &Vec<String>, labels: &Vec<Span>) -> Vec<usize> {
    let mut motif_counts = vec![0; motifs.len()];
    for span in labels {
        motif_counts[span.motif_index] += 1;
    }
    motif_counts
}

pub fn collapse_labels(spans: Vec<Span>) -> Vec<Span> {
    let mut collapsed = Vec::new();
    for span in spans {
        if collapsed.is_empty() {
            collapsed.push(span);
            continue;
        }

        let last_span = collapsed.last_mut().unwrap();
        if last_span.motif_index == span.motif_index && last_span.end == span.start {
            last_span.end = span.end;
        } else {
            collapsed.push(span);
        }
    }
    collapsed
}

pub fn replace_invalid_bases(seq: &str, allowed_bases: &[char]) -> String {
    seq.as_bytes()
        .iter()
        .enumerate()
        .map(|(index, base)| {
            let base = *base as char;
            if allowed_bases.contains(&base) {
                base
            } else {
                allowed_bases[index % allowed_bases.len()]
            }
        })
        .collect()
}
