use super::events::{get_events, HmmEvent};
use super::hmm::Hmm;

pub fn calc_purity(query: &[u8], hmm: &Hmm, motifs: &[Vec<u8>], states: &[usize]) -> f64 {
    if query.is_empty() {
        return f64::NAN;
    }

    let events = get_events(hmm, motifs, states, query);
    let edit_dist = events
        .iter()
        .filter(|e| {
            [
                HmmEvent::Del,
                HmmEvent::Ins,
                HmmEvent::Mismatch,
                HmmEvent::Skip,
            ]
            .contains(*e)
        })
        .count() as f64;

    // Length of the implicit reference against which the query is aligned
    let ref_len = events
        .iter()
        .filter(|e| {
            [
                HmmEvent::Match,
                HmmEvent::Mismatch,
                HmmEvent::Del,
                HmmEvent::Skip,
            ]
            .contains(*e)
        })
        .count();

    let max_dist = std::cmp::max(ref_len, query.len()) as f64;
    (max_dist - edit_dist) / max_dist
}

#[cfg(test)]
mod tests {
    use super::super::builder::build_hmm;
    use super::*;

    #[test]
    fn calculate_purity_of_perfect_repeats() {
        let motifs = vec!["CAG".as_bytes().to_vec(), "CCG".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        let query = "CAGCAGCAGCCGCCGCCGCCG";
        let states = hmm.label(&query);
        assert_eq!(calc_purity(&query.as_bytes(), &hmm, &motifs, &states), 1.0);
    }

    #[test]
    fn calculate_purity_of_imperfect_repeats() {
        let motifs = vec!["CAG".as_bytes().to_vec(), "CCG".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        let query = "CAGCGCAGCCGCCGCCGGG";
        let states = hmm.label(&query);

        let purity = calc_purity(&query.as_bytes(), &hmm, &motifs, &states);
        assert_eq!(purity, 18.0 / 21.0);
    }

    #[test]
    fn calculate_purity_of_repeats_with_skip_states() {
        let motifs = vec!["CAG".as_bytes().to_vec(), "CCG".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        let query = "CAGCAGCAGTTTTTTTTCCGCCGCCG";
        let states = hmm.label(&query);

        let purity = calc_purity(&query.as_bytes(), &hmm, &motifs, &states);
        assert_eq!(purity, 18.0 / 26.0);
    }

    #[test]
    fn calculate_purity_of_empty_query() {
        let motifs = vec!["CAG".as_bytes().to_vec(), "CCG".as_bytes().to_vec()];
        let hmm = build_hmm(&motifs);
        let query = "";
        let states = hmm.label(&query);
        let purity = calc_purity(&query.as_bytes(), &hmm, &motifs, &states);
        assert!(purity.is_nan());
    }
}
