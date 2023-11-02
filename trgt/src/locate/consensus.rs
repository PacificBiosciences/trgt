use arrayvec::ArrayVec;

pub fn get_consensus(sizes: ArrayVec<usize, 2>, seqs: &[&str], counts: &[usize]) -> Vec<String> {
    let mut consensuses = Vec::new();

    let allele = get_closest_size(seqs, sizes[0]).unwrap();
    let consensus = get_most_frequent_seq(seqs, counts, allele);
    consensuses.push(consensus);

    if sizes.len() != 1 && sizes[0] != sizes[1] {
        let allele = get_closest_size(seqs, sizes[1]).unwrap();
        let consensus = get_most_frequent_seq(seqs, counts, allele);
        consensuses.push(consensus);
    }

    consensuses
}

fn get_closest_size(seqs: &[&str], allele: usize) -> Option<usize> {
    let mut closest_size = None;
    for seq in seqs {
        let read_len = seq.len();

        if closest_size.is_none() {
            closest_size = Some(read_len);
            continue;
        }

        if closest_size.unwrap().abs_diff(allele) > read_len.abs_diff(allele) {
            closest_size = Some(read_len);
        }
    }

    closest_size
}

fn get_most_frequent_seq(seqs: &[&str], counts: &[usize], length: usize) -> String {
    seqs.iter()
        .zip(counts)
        .filter(|rec| rec.0.len() == length)
        .max_by_key(|rec| rec.1)
        .unwrap()
        .0
        .to_string()

    //   let mut hash: HashMap<&str, usize> = HashMap::new();
    //
    // for seq in seqs {
    //     if seq.len() == length {
    //         *hash.entry(*seq).or_insert(0) += 1;
    //      }
    //   }
    //
    //let consensus = hash.iter().max_by_key(|rec| rec.1).unwrap().0;
    //consensus.to_string()
}
