#[derive(Debug, PartialEq, Clone)]
pub enum Mult {
    Once,
    Many,
}

#[derive(Debug, PartialEq, Clone)]
pub struct Seq {
    pub motif: String,
    pub mult: Mult,
}

pub type Struc = Vec<Seq>;

#[derive(Debug, PartialEq, Clone)]
pub enum RegionLabel {
    Flank(usize, usize),      // Coordinates
    Tr(usize, usize, String), // Coordinates, Motif
    Seq(usize, usize),
    Other(usize, usize),
}

pub fn struc(encoding: &str) -> Struc {
    let mut struc = Vec::new();
    for piece in encoding.split(|c| c == '(' || c == 'n') {
        if piece.is_empty() {
            continue;
        }

        if piece.len() >= 2 && &piece[piece.len() - 1..] == ")" {
            let seq = piece[..piece.len() - 1].to_string();
            let mult = Mult::Many;
            struc.push(Seq { motif: seq, mult });
        } else {
            struc.push(Seq {
                motif: piece.to_string(),
                mult: Mult::Once,
            });
        }
    }

    struc
}

#[cfg(test)]
mod tests {
    use crate::struc::*;

    #[test]
    fn test_decoding_homopolymer_struc() {
        let struc = struc("(A)n");

        assert_eq!(
            struc,
            vec![Seq {
                motif: "A".to_string(),
                mult: Mult::Many,
            }]
        );
    }

    #[test]
    fn test_decoding_simple_struc() {
        let struc = struc("(CAG)n");

        assert_eq!(
            struc,
            vec![Seq {
                motif: "CAG".to_string(),
                mult: Mult::Many,
            }]
        );
    }

    #[test]
    fn test_decoding_complex_struc() {
        let struc = struc("(CAG)nCAACAG(CCG)n");

        assert_eq!(
            struc,
            vec![
                Seq {
                    motif: "CAG".to_string(),
                    mult: Mult::Many,
                },
                Seq {
                    motif: "CAACAG".to_string(),
                    mult: Mult::Once,
                },
                Seq {
                    motif: "CCG".to_string(),
                    mult: Mult::Many,
                }
            ]
        );
    }
}
