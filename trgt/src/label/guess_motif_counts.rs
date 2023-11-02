use crate::label::struc::*;

#[derive(Debug, Clone)]
struct Node {
    start: usize,
    end: usize,
    span: usize,
    pred: Option<(usize, usize)>,
    block_index: usize,
}

type MotifGraph = Vec<Vec<Node>>;

pub fn guess_motif_counts(motifs: &Vec<Motif>, query: &str) -> Vec<usize> {
    let graph = make_graph(motifs, query);
    count_motifs(motifs, &graph)
}

fn make_graph(motifs: &Vec<Motif>, query: &str) -> MotifGraph {
    let mut graph = Vec::with_capacity(motifs.len());
    graph.resize(motifs.len(), vec![]);

    for kmer_start in 0..query.len() {
        for (seq_index, seq) in motifs.iter().enumerate() {
            let motif = &seq.seq;
            let kmer_end = kmer_start + motif.len();
            if kmer_end > query.len() {
                continue;
            }
            let kmer = &query[kmer_start..kmer_end];
            if kmer == motif {
                let pred = get_pred(&graph, seq_index, kmer_start);
                let span = match pred {
                    None => kmer.len(),
                    Some((seq_index, node_index)) => kmer.len() + graph[seq_index][node_index].span,
                };

                graph[seq_index].push(Node {
                    start: kmer_start,
                    end: kmer_end,
                    span,
                    pred,
                    block_index: seq_index,
                });
            }
        }
    }

    graph
}

fn get_pred(
    graph: &MotifGraph,
    target_seq_index: usize,
    target_start: usize,
) -> Option<(usize, usize)> {
    let mut pred = None;

    for seq_index in (0..=target_seq_index).rev() {
        let seq_nodes = &graph[seq_index];
        for node_index in (0..seq_nodes.len()).rev() {
            let node = &seq_nodes[node_index];
            if node.end <= target_start {
                pred = match pred {
                    None => Some((seq_index, node_index)),
                    Some((prev_seq_index, prev_node_index)) => {
                        let prev_node = &graph[prev_seq_index][prev_node_index];
                        if node.span > prev_node.span {
                            Some((seq_index, node_index))
                        } else {
                            pred
                        }
                    }
                };
                break;
            }
        }
    }

    pred
}

fn count_motifs(motifs: &Vec<Motif>, graph: &MotifGraph) -> Vec<usize> {
    // Find traceback index
    let mut last_indexes = None;

    for (motif_index, motif_nodes) in graph.iter().enumerate() {
        if motif_nodes.is_empty() {
            continue;
        }
        let node_index = motif_nodes
            .iter()
            .enumerate()
            .max_by(|x, y| x.1.span.cmp(&y.1.span))
            .unwrap()
            .0;

        let candidate_motif = &graph[motif_index][node_index];

        last_indexes = match last_indexes {
            None => Some((motif_index, node_index)),
            Some((last_seq_index, last_node_index)) => {
                let last_block = &graph[last_seq_index][last_node_index];
                if candidate_motif.span > last_block.span {
                    Some((motif_index, node_index))
                } else {
                    last_indexes
                }
            }
        }
    }

    // Perform traceback
    match last_indexes {
        Some(indexes) => traceback(motifs, graph, indexes),
        None => {
            vec![0; motifs.len()]
        }
    }
}

fn traceback(motifs: &Vec<Motif>, graph: &MotifGraph, last_indexes: (usize, usize)) -> Vec<usize> {
    let mut motif_spans = Vec::new();
    motif_spans.resize(motifs.len(), None);
    let (last_block, last_node) = last_indexes;
    let mut node = &graph[last_block][last_node];

    loop {
        if motif_spans[node.block_index].is_none() {
            motif_spans[node.block_index] = Some((node.start, node.end));
        } else {
            let (start, end) = motif_spans[node.block_index].unwrap();
            assert!(node.end <= start);
            motif_spans[node.block_index] = Some((node.start, end));
        }

        if node.pred.is_none() {
            break;
        } else {
            node = &graph[node.pred.unwrap().0][node.pred.unwrap().1];
        }
    }

    let mut motif_counts = Vec::new();
    for (index, span) in motif_spans.iter().enumerate() {
        let motif = &motifs[index].seq;
        motif_counts.push(match span {
            None => 0,
            Some((start, end)) => ((end - start) as f64 / motif.len() as f64).round() as usize,
        });
    }

    motif_counts
}

#[cfg(test)]
mod tests {
    use crate::label::guess_motif_counts::*;

    #[test]
    fn test_matching_simple_patterns_with_kmers() {
        let motifs = decode_regexp("(CCG)n");

        assert_eq!(guess_motif_counts(&motifs, "CCGCCGCCGCCG"), vec![4]);
        assert_eq!(guess_motif_counts(&motifs, "CCGCCGCGCCG"), vec![4]);
        assert_eq!(guess_motif_counts(&motifs, "ATATATATATA"), vec![0]);
    }

    #[test]
    fn test_matching_complex_patterns_with_kmers() {
        let motifs = decode_regexp("(CAG)nCAACAG(CCG)n");

        let query = "CAGCAGCAGCAACAGCCGCCG";
        assert_eq!(guess_motif_counts(&motifs, query), vec![3, 1, 2]);
    }

    #[test]
    fn test_matching_htt_repeat_with_missing_interruption() {
        let motifs = decode_regexp("(CAG)nCAACAG(CCG)n");

        //                 111111111111111333333
        let query = "CAGCAGCAGCATCAGCCGCCG";
        assert_eq!(guess_motif_counts(&motifs, query), vec![5, 0, 2]);
    }
}
