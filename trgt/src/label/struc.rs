#[derive(Debug, PartialEq)]
pub enum Mult {
    Once,
    Many,
}

#[derive(Debug, PartialEq)]
pub struct Motif {
    pub seq: String,
    pub mult: Mult,
}

#[derive(Debug, PartialEq)]
pub enum TrType {
    RegExp,
    SingleMotif,
    Hmm,
}

#[derive(Debug, PartialEq)]
pub struct TrDef {
    pub motifs: Vec<Motif>,
    pub tr_type: TrType,
}

pub fn decode_regexp(encoding: &str) -> Vec<Motif> {
    let mut motifs = Vec::new();
    for piece in encoding.split(|c| c == '(' || c == 'n') {
        if piece.is_empty() {
            continue;
        }

        if piece.len() >= 2 && &piece[piece.len() - 1..] == ")" {
            let seq = piece[..piece.len() - 1].to_string();
            let mult = Mult::Many;
            motifs.push(Motif { seq, mult });
        } else {
            motifs.push(Motif {
                seq: piece.to_string(),
                mult: Mult::Once,
            });
        }
    }

    motifs
}

/*
pub fn decode_motif_list(encoding: &str) -> Vec<Motif> {
    let mut motifs = Vec::new();
    for piece in encoding.split(",") {
        if piece.is_empty() {
            continue;
        }
        motifs.push(Motif {
            seq: piece.to_string(),
            mult: Mult::Many,
        });
    }
    motifs
} */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cag_repeat() {
        let motifs = decode_regexp("(CAG)n");

        let expected = vec![Motif {
            seq: "CAG".to_string(),
            mult: Mult::Many,
        }];
        assert_eq!(motifs, expected);
    }

    #[test]
    fn test_htt_repeat() {
        let motifs = decode_regexp("(CAG)nCAACAG(CCG)n");

        let expected = vec![
            Motif {
                seq: "CAG".to_string(),
                mult: Mult::Many,
            },
            Motif {
                seq: "CAACAG".to_string(),
                mult: Mult::Once,
            },
            Motif {
                seq: "CCG".to_string(),
                mult: Mult::Many,
            },
        ];
        assert_eq!(motifs, expected);
    }

    #[test]
    fn test_fxn_repeat() {
        let motifs = decode_regexp("(A)n(GAA)n");

        let expected = vec![
            Motif {
                seq: "A".to_string(),
                mult: Mult::Many,
            },
            Motif {
                seq: "GAA".to_string(),
                mult: Mult::Many,
            },
        ];
        assert_eq!(motifs, expected);
    }
}
