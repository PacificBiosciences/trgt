// https://github.com/rust-bio/rust-htslib/blob/master/src/faidx/mod.rs
// Copyright 2020 Manuel Landesfeind, Evotec International GmbH
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

// TODO: Temporary work-around until https://github.com/rust-bio/rust-htslib/pull/410 gets merged which adds fetch_seq_len

//!
//! Module for working with faidx-indexed FASTA files.
//!

use rust_htslib::errors::{Error, Result};
use rust_htslib::htslib;
use rust_htslib::utils::path_as_bytes;
use std::collections::HashMap;
use std::ffi;
use std::path::Path;

/// A Fasta reader.
#[derive(Debug)]
pub struct Reader {
    inner: *mut htslib::faidx_t,
}

impl Reader {
    /// Create a new Reader from a path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, Error> {
        Self::new(&path_as_bytes(path, true)?)
    }

    /// Internal function to create a Reader from some sort of path (could be file path but also URL).
    /// The path or URL will be handled by the c-implementation transparently.
    ///
    /// # Arguments
    ///
    /// * `path` - the path or URL to open
    fn new(path: &[u8]) -> Result<Self, Error> {
        let cpath = ffi::CString::new(path).unwrap();
        let inner = unsafe { htslib::fai_load(cpath.as_ptr()) };
        Ok(Self { inner })
    }

    /// Fetch the sequence as a byte array.
    ///
    /// # Arguments
    ///
    /// * `name` - the name of the template sequence (e.g., "chr1")
    /// * `begin` - the offset within the template sequence (starting with 0)
    /// * `end` - the end position to return (if smaller than `begin`, the behavior is undefined).
    pub fn fetch_seq<N: AsRef<str>>(&self, name: N, begin: usize, end: usize) -> Result<&[u8]> {
        if begin > std::i64::MAX as usize {
            return Err(Error::FaidxPositionTooLarge);
        }
        if end > std::i64::MAX as usize {
            return Err(Error::FaidxPositionTooLarge);
        }
        let cname = ffi::CString::new(name.as_ref().as_bytes()).unwrap();
        let len_out: i64 = 0;
        let cseq = unsafe {
            let ptr = htslib::faidx_fetch_seq64(
                self.inner,                          //*const faidx_t,
                cname.as_ptr(),                      // c_name
                begin as htslib::hts_pos_t,          // p_beg_i
                end as htslib::hts_pos_t,            // p_end_i
                &mut (len_out as htslib::hts_pos_t), //len
            );
            ffi::CStr::from_ptr(ptr)
        };

        Ok(cseq.to_bytes())
    }

    /// Fetches the sequence and returns it as string.
    ///
    /// # Arguments
    ///
    /// * `name` - the name of the template sequence (e.g., "chr1")
    /// * `begin` - the offset within the template sequence (starting with 0)
    /// * `end` - the end position to return (if smaller than `begin`, the behavior is undefined).
    pub fn fetch_seq_string<N: AsRef<str>>(
        &self,
        name: N,
        begin: usize,
        end: usize,
    ) -> Result<String> {
        let bytes = self.fetch_seq(name, begin, end)?;
        Ok(std::str::from_utf8(bytes).unwrap().to_owned())
    }

    /// Fetches the number of sequences in the fai index
    pub fn n_seqs(&self) -> u64 {
        let n = unsafe { htslib::faidx_nseq(self.inner) };
        n as u64
    }

    /// Fetches the i-th sequence name
    ///
    /// # Arguments
    ///
    /// * `i` - index to query
    pub fn seq_name(&self, i: i32) -> Result<String> {
        let cname = unsafe {
            let ptr = htslib::faidx_iseq(self.inner, i);
            ffi::CStr::from_ptr(ptr)
        };

        let out = match cname.to_str() {
            Ok(s) => s.to_string(),
            Err(_) => {
                return Err(Error::FaidxBadSeqName);
            }
        };

        Ok(out)
    }

    pub fn fetch_seq_len<N: AsRef<str>>(&self, name: N) -> Option<i32> {
        let cname = ffi::CString::new(name.as_ref().as_bytes()).unwrap();
        let seq_len = unsafe { htslib::faidx_seq_len(self.inner, cname.as_ptr()) };
        if seq_len >= 0 {
            Some(seq_len)
        } else {
            None
        }
    }

    /// Create a HashMap mapping each sequence name to its length.
    pub fn create_chrom_lookup(&self) -> Result<HashMap<String, u32>, String> {
        let num_seqs = self.n_seqs() as usize;
        let mut map = HashMap::with_capacity(num_seqs);
        for i in 0..num_seqs {
            let name = self.seq_name(i as i32).map_err(|e| e.to_string())?;
            if let Some(len) = self.fetch_seq_len(&name) {
                let len_u32 = u32::try_from(len).map_err(|_| {
                    format!(
                        "Sequence length for '{}' is negative and cannot be converted to u32",
                        &name
                    )
                })?;
                map.insert(name, len_u32);
            }
        }
        Ok(map)
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        unsafe {
            htslib::fai_destroy(self.inner);
        }
    }
}
