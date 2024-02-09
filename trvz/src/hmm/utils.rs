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

pub fn replace_invalid_bases(seq: &str) -> String {
    seq.as_bytes()
        .iter()
        .enumerate()
        .map(|(i, b)| match b {
            b'A' | b'T' | b'C' | b'G' => *b as char,
            _ => ['A', 'T', 'C', 'G'][i % 4],
        })
        .collect()
}
