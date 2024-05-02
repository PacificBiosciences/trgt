use super::read::MethInfo;

pub fn decode_on_plus(bases: &[u8], meth: &MethInfo) -> Option<Vec<u8>> {
    let mut profile = Vec::new();

    let mut num_cs_skipped = 0;
    let mut cite_index = 0;
    let mut non_cpgs_called = 0;

    for (index, base) in bases.iter().enumerate() {
        if *base != b'C' || index + 1 == bases.len() {
            continue;
        }

        let dinuc = &bases[index..index + 2];

        if dinuc == b"CG" {
            profile.push(0);
        }

        // TODO: Check if the first condition needed
        if cite_index != meth.poses.len() && num_cs_skipped == meth.poses[cite_index] {
            if dinuc == b"CG" {
                *profile.last_mut().unwrap() = meth.probs[cite_index];
            } else {
                non_cpgs_called += 1;
            }

            num_cs_skipped = 0;
            cite_index += 1;
        } else {
            num_cs_skipped += 1;
        }
    }

    if non_cpgs_called > 0 {
        log::warn!("Warning: non_cpgs_called = {non_cpgs_called}");
    }

    Some(profile)
}

pub fn decode_on_minus(bases: &[u8], meth: &MethInfo) -> Option<Vec<u8>> {
    let mut profile = Vec::new();

    let mut num_cs_skipped = 0;
    let mut cite_index = 0;
    let mut non_cpgs_called = 0;

    for (index_rev, base) in bases.iter().rev().enumerate() {
        let index = bases.len() - 1 - index_rev;

        if *base != b'G' || index == 0 {
            continue;
        }

        let dinuc = &bases[index - 1..index + 1];

        if dinuc == b"CG" {
            profile.push(0);
        }

        if cite_index != meth.poses.len() && num_cs_skipped == meth.poses[cite_index] {
            if dinuc == b"CG" {
                *profile.last_mut().unwrap() = meth.probs[cite_index];
            } else {
                non_cpgs_called += 1;
            }

            num_cs_skipped = 0;
            cite_index += 1;
        } else {
            num_cs_skipped += 1;
        }
    }

    if non_cpgs_called > 0 {
        log::warn!("Warning: non_cpgs_called = {non_cpgs_called}");
    }

    Some(profile)
}
