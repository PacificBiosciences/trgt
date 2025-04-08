use crate::wfa2;
use bio::alignment::AlignmentOperation;
use core::slice;
use std::ptr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryModel {
    MemoryHigh,
    MemoryMed,
    MemoryLow,
    MemoryUltraLow,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentScope {
    Score,
    Alignment,
}

impl From<u32> for AlignmentScope {
    fn from(value: u32) -> Self {
        match value {
            wfa2::alignment_scope_t_compute_score => AlignmentScope::Score,
            wfa2::alignment_scope_t_compute_alignment => AlignmentScope::Alignment,
            _ => panic!("Unknown alignment scope: {}", value),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Heuristic {
    None,
    WFadaptive(i32, i32, i32),
    WFmash(i32, i32, i32),
    XDrop(i32, i32),
    ZDrop(i32, i32),
    BandedStatic(i32, i32),
    BandedAdaptive(i32, i32, i32),
}

pub enum DistanceMetric {
    Indel,
    Edit,
    GapLinear,
    GapAffine,
    GapAffine2p,
}

impl From<u32> for DistanceMetric {
    fn from(value: u32) -> Self {
        match value {
            wfa2::distance_metric_t_indel => DistanceMetric::Indel,
            wfa2::distance_metric_t_edit => DistanceMetric::Edit,
            wfa2::distance_metric_t_gap_linear => DistanceMetric::GapLinear,
            wfa2::distance_metric_t_gap_affine => DistanceMetric::GapAffine,
            wfa2::distance_metric_t_gap_affine_2p => DistanceMetric::GapAffine2p,
            _ => panic!("Unknown distance metric: {}", value),
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum AlignmentStatus {
    // OK Status (>=0)
    StatusAlgCompleted = wfa2::WF_STATUS_ALG_COMPLETED as isize,
    StatusAlgPartial = wfa2::WF_STATUS_ALG_PARTIAL as isize,
    // FAILED Status (<0)
    StatusMaxStepsReached = wfa2::WF_STATUS_MAX_STEPS_REACHED as isize,
    StatusOOM = wfa2::WF_STATUS_OOM as isize,
    StatusUnattainable = wfa2::WF_STATUS_UNATTAINABLE as isize,
}

impl From<i32> for AlignmentStatus {
    fn from(value: i32) -> Self {
        match value {
            x if x == wfa2::WF_STATUS_ALG_COMPLETED as i32 => AlignmentStatus::StatusAlgCompleted,
            x if x == wfa2::WF_STATUS_ALG_PARTIAL as i32 => AlignmentStatus::StatusAlgPartial,
            wfa2::WF_STATUS_MAX_STEPS_REACHED => AlignmentStatus::StatusMaxStepsReached,
            wfa2::WF_STATUS_OOM => AlignmentStatus::StatusOOM,
            wfa2::WF_STATUS_UNATTAINABLE => AlignmentStatus::StatusUnattainable,
            _ => panic!("Unknown alignment status: {}", value),
        }
    }
}

#[derive(Debug, Copy, Clone)]
struct WFAttributes {
    inner: wfa2::wavefront_aligner_attr_t,
}

impl WFAttributes {
    fn default() -> Self {
        Self {
            inner: unsafe { wfa2::wavefront_aligner_attr_default },
        }
    }

    fn memory_model(mut self, memory_model: MemoryModel) -> Self {
        let memory_mode = match memory_model {
            MemoryModel::MemoryHigh => wfa2::wavefront_memory_t_wavefront_memory_high,
            MemoryModel::MemoryMed => wfa2::wavefront_memory_t_wavefront_memory_med,
            MemoryModel::MemoryLow => wfa2::wavefront_memory_t_wavefront_memory_low,
            MemoryModel::MemoryUltraLow => wfa2::wavefront_memory_t_wavefront_memory_ultralow,
        };
        self.inner.memory_mode = memory_mode;
        self
    }

    fn alignment_scope(mut self, alignment_scope: AlignmentScope) -> Self {
        let alignment_scope = match alignment_scope {
            AlignmentScope::Score => wfa2::alignment_scope_t_compute_score,
            AlignmentScope::Alignment => wfa2::alignment_scope_t_compute_alignment,
        };
        self.inner.alignment_scope = alignment_scope;
        self
    }

    fn indel_penalties(mut self) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_indel;
        self
    }

    fn edit_penalties(mut self) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_edit;
        self
    }

    fn linear_penalties(mut self, match_: i32, mismatch: i32, indel: i32) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_gap_linear;
        self.inner.linear_penalties.match_ = match_; // (Penalty representation usually M <= 0)
        self.inner.linear_penalties.mismatch = mismatch; // (Penalty representation usually X > 0)
        self.inner.linear_penalties.indel = indel; // (Penalty representation usually I > 0)
        self
    }

    fn affine_penalties(
        mut self,
        match_: i32,
        mismatch: i32,
        gap_opening: i32,
        gap_extension: i32,
    ) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_gap_affine;
        self.inner.affine_penalties.match_ = match_; // (Penalty representation usually M <= 0)
        self.inner.affine_penalties.mismatch = mismatch; // (Penalty representation usually X > 0)
        self.inner.affine_penalties.gap_opening = gap_opening; // (Penalty representation usually O > 0)
        self.inner.affine_penalties.gap_extension = gap_extension; // (Penalty representation usually E > 0)
        self
    }

    fn affine2p_penalties(
        mut self,
        match_: i32,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
    ) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_gap_affine_2p;
        self.inner.affine2p_penalties.match_ = match_; // (Penalty representation usually M <= 0)
        self.inner.affine2p_penalties.mismatch = mismatch; // (Penalty representation usually X > 0)
                                                           // Usually concave Q1 + E1 < Q2 + E2 and E1 > E2.
        self.inner.affine2p_penalties.gap_opening1 = gap_opening1; // (Penalty representation usually O1 > 0)
        self.inner.affine2p_penalties.gap_extension1 = gap_extension1; // (Penalty representation usually E1 > 0)
        self.inner.affine2p_penalties.gap_opening2 = gap_opening2; // (Penalty representation usually O2 > 0)
        self.inner.affine2p_penalties.gap_extension2 = gap_extension2; // (Penalty representation usually E2 > 0)
        self
    }
}

pub struct WFAligner {
    attributes: WFAttributes,
    inner: *mut wfa2::wavefront_aligner_t,
}

impl WFAligner {
    pub fn new(alignment_scope: AlignmentScope, memory_model: MemoryModel) -> Self {
        let attributes = WFAttributes::default()
            .memory_model(memory_model)
            .alignment_scope(alignment_scope);
        Self {
            attributes,
            inner: ptr::null_mut(),
        }
    }

    /// Returns the penalties as a tuple based on the distance metric.
    pub fn get_penalties(&self) -> (i32, i32, Option<i32>, Option<i32>, Option<i32>, Option<i32>) {
        match DistanceMetric::from(self.attributes.inner.distance_metric) {
            DistanceMetric::Edit => (
                self.attributes.inner.linear_penalties.match_,
                self.attributes.inner.linear_penalties.mismatch,
                Some(self.attributes.inner.linear_penalties.indel),
                None,
                None,
                None,
            ),
            DistanceMetric::Indel => (
                self.attributes.inner.linear_penalties.match_,
                self.attributes.inner.linear_penalties.mismatch,
                Some(self.attributes.inner.linear_penalties.indel),
                None,
                None,
                None,
            ),
            DistanceMetric::GapLinear => (
                self.attributes.inner.linear_penalties.match_,
                self.attributes.inner.linear_penalties.mismatch,
                Some(self.attributes.inner.linear_penalties.indel),
                None,
                None,
                None,
            ),
            DistanceMetric::GapAffine => (
                self.attributes.inner.affine_penalties.match_,
                self.attributes.inner.affine_penalties.mismatch,
                Some(self.attributes.inner.affine_penalties.gap_opening),
                Some(self.attributes.inner.affine_penalties.gap_extension),
                None,
                None,
            ),
            DistanceMetric::GapAffine2p => (
                self.attributes.inner.affine2p_penalties.match_,
                self.attributes.inner.affine2p_penalties.mismatch,
                Some(self.attributes.inner.affine2p_penalties.gap_opening1),
                Some(self.attributes.inner.affine2p_penalties.gap_extension1),
                Some(self.attributes.inner.affine2p_penalties.gap_opening2),
                Some(self.attributes.inner.affine2p_penalties.gap_extension2),
            ),
        }
    }
}

impl Drop for WFAligner {
    fn drop(&mut self) {
        unsafe {
            if !self.inner.is_null() {
                wfa2::wavefront_aligner_delete(self.inner);
            }
        }
    }
}

impl WFAligner {
    pub fn set_alignment_end_to_end(&mut self) {
        unsafe {
            wfa2::wavefront_aligner_set_alignment_end_to_end(self.inner);
        }
    }

    pub fn set_alignment_ends_free(
        &mut self,
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text_begin_free: i32,
        text_end_free: i32,
    ) {
        unsafe {
            wfa2::wavefront_aligner_set_alignment_free_ends(
                self.inner,
                pattern_begin_free,
                pattern_end_free,
                text_begin_free,
                text_end_free,
            );
        }
    }

    pub fn align_end_to_end(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentStatus {
        let status = unsafe {
            self.set_alignment_end_to_end();
            wfa2::wavefront_align(
                self.inner,
                pattern.as_ptr() as *const i8,
                pattern.len() as i32,
                text.as_ptr() as *const i8,
                text.len() as i32,
            )
        };
        AlignmentStatus::from(status)
    }

    pub fn align_ends_free(
        &mut self,
        pattern: &[u8],
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text: &[u8],
        text_begin_free: i32,
        text_end_free: i32,
    ) -> AlignmentStatus {
        let status = unsafe {
            self.set_alignment_ends_free(
                pattern_begin_free,
                pattern_end_free,
                text_begin_free,
                text_end_free,
            );
            wfa2::wavefront_align(
                self.inner,
                pattern.as_ptr() as *const i8,
                pattern.len() as i32,
                text.as_ptr() as *const i8,
                text.len() as i32,
            )
        };
        AlignmentStatus::from(status)
    }

    pub fn score(&self) -> i32 {
        unsafe { *(*self.inner).cigar }.score
    }

    fn cigar_score_gap_affine2p_get_operations_score(
        operation: char,
        op_length: i32,
        penalties: &wfa2::affine2p_penalties_t,
    ) -> i32 {
        match operation {
            'M' => op_length * penalties.match_,
            'X' => op_length * penalties.mismatch,
            'D' | 'I' => {
                let score1 = penalties.gap_opening1 + penalties.gap_extension1 * op_length;
                let score2 = penalties.gap_opening2 + penalties.gap_extension2 * op_length;
                std::cmp::min(score1, score2)
            }
            _ => panic!("Invalid operation: {}", operation),
        }
    }

    fn cigar_score_gap_affine2_clipped(&self, flank_len: usize) -> i32 {
        let cigar = unsafe { (*self.inner).cigar.as_ref() }.unwrap();
        let begin_offset = cigar.begin_offset as isize + flank_len as isize;
        let end_offset = cigar.end_offset as isize - flank_len as isize;

        let operations: Vec<char> =
            unsafe { std::slice::from_raw_parts(cigar.operations, cigar.max_operations as usize) }
                .iter()
                .map(|&op| op as u8 as char)
                .collect();

        let mut score = 0;
        let mut op_length = 0;
        let mut last_op: Option<char> = None;
        for i in begin_offset..end_offset {
            let cur_op = operations[i as usize];
            if let Some(op) = last_op {
                if cur_op != op {
                    score -= Self::cigar_score_gap_affine2p_get_operations_score(
                        op,
                        op_length,
                        &self.attributes.inner.affine2p_penalties,
                    );
                    op_length = 0;
                }
            }
            last_op = Some(cur_op);
            op_length += 1;
        }
        if let Some(op) = last_op {
            score -= Self::cigar_score_gap_affine2p_get_operations_score(
                op,
                op_length,
                &self.attributes.inner.affine2p_penalties,
            );
        }
        score
    }

    pub fn cigar_score_clipped(&self, flank_len: usize) -> i32 {
        if AlignmentScope::from(self.attributes.inner.alignment_scope) == AlignmentScope::Score {
            panic!("Cannot clip when AlignmentScope is Score");
        }

        match DistanceMetric::from(self.attributes.inner.distance_metric) {
            DistanceMetric::Indel | DistanceMetric::Edit => {
                unimplemented!("Indel/Linear distance metric not implemented")
            }
            DistanceMetric::GapLinear => {
                unimplemented!("GapLinear distance metric not implemented")
            }
            DistanceMetric::GapAffine => {
                unimplemented!("GapAffine distance metric not implemented")
            }
            DistanceMetric::GapAffine2p => self.cigar_score_gap_affine2_clipped(flank_len),
        }
    }

    pub fn get_sam_cigar(&self, show_mismatches: bool) -> Vec<u32> {
        if AlignmentScope::from(self.attributes.inner.alignment_scope) == AlignmentScope::Score {
            panic!("Cannot get SAM CIGAR when AlignmentScope is Score");
        }
        if self.inner.is_null() {
            panic!("Internal aligner pointer is null");
        }

        unsafe {
            let mut sam_cigar_buffer_ptr: *mut u32 = std::ptr::null_mut();
            let mut sam_cigar_length: i32 = 0;

            wfa2::cigar_get_CIGAR(
                (*self.inner).cigar,
                show_mismatches,
                &mut sam_cigar_buffer_ptr,
                &mut sam_cigar_length,
            );

            if !sam_cigar_buffer_ptr.is_null() && sam_cigar_length > 0 {
                let cigar_buffer_slice =
                    slice::from_raw_parts(sam_cigar_buffer_ptr, sam_cigar_length as usize);
                cigar_buffer_slice.to_vec()
            } else {
                Vec::new()
            }
        }
    }

    /// Decodes a raw SAM CIGAR buffer (Vec<u32>) into a vector of (length, operation) tuples.
    pub fn decode_sam_cigar(sam_cigar_buffer: &[u32]) -> Vec<(u32, char)> {
        sam_cigar_buffer
            .iter()
            .map(|&encoded_op| {
                let len = encoded_op >> 4; // Length is in the upper 28 bits
                let op_code = encoded_op & 0xF; // Operation code is in the lower 4 bits
                let op_char = match op_code {
                    0 => 'M', // BAM_CMATCH (Alignment match (can be sequence match or mismatch))
                    1 => 'I', // BAM_CINS (Insertion to the reference)
                    2 => 'D', // BAM_CDEL (Deletion from the reference)
                    3 => 'N', // BAM_CREF_SKIP (Skipped region from the reference)
                    4 => 'S', // BAM_CSOFT_CLIP (Soft clipping (clipped sequences present in SEQ))
                    5 => 'H', // BAM_CHARD_CLIP (Hard clipping (clipped sequences NOT present in SEQ))
                    6 => 'P', // BAM_CPAD (Padding (silent deletion from padded reference))
                    7 => '=', // BAM_CEQUAL (Sequence match)
                    8 => 'X', // BAM_CDIFF (Sequence mismatch)
                    _ => '?', // Unknown operation
                };
                (len, op_char)
            })
            .collect()
    }

    /// Counts the number of match ('M') operations in the CIGAR string.
    pub fn count_matches(&self) -> i32 {
        if AlignmentScope::from(self.attributes.inner.alignment_scope) == AlignmentScope::Score {
            panic!("Cannot count matches when AlignmentScope is Score");
        }
        if self.inner.is_null() {
            panic!("Internal aligner pointer is null");
        }
        let cigar_ptr = unsafe { (*self.inner).cigar };
        if cigar_ptr.is_null() {
            panic!("CIGAR pointer is null, cannot count matches.");
        }
        unsafe { wfa2::cigar_count_matches(cigar_ptr) }
    }

    pub fn cigar_score(&mut self) -> i32 {
        if AlignmentScope::from(self.attributes.inner.alignment_scope) == AlignmentScope::Score {
            panic!("Cannot calculate CIGAR score when AlignmentScope is Score");
        }

        unsafe {
            match DistanceMetric::from(self.attributes.inner.distance_metric) {
                DistanceMetric::Indel | DistanceMetric::Edit => {
                    wfa2::cigar_score_edit((*self.inner).cigar)
                }
                DistanceMetric::GapLinear => wfa2::cigar_score_gap_linear(
                    (*self.inner).cigar,
                    &mut self.attributes.inner.linear_penalties,
                ),
                DistanceMetric::GapAffine => wfa2::cigar_score_gap_affine(
                    (*self.inner).cigar,
                    &mut self.attributes.inner.affine_penalties,
                ),
                DistanceMetric::GapAffine2p => wfa2::cigar_score_gap_affine2p(
                    (*self.inner).cigar,
                    &mut self.attributes.inner.affine2p_penalties,
                ),
            }
        }
    }

    pub fn cigar(&self, flank_len: Option<usize>) -> String {
        let offset = flank_len.unwrap_or(0);
        let mut cstr = String::new();

        let cigar = unsafe { (*self.inner).cigar.as_ref() }.unwrap();

        let begin_offset = cigar.begin_offset as usize + offset;
        let end_offset = cigar.end_offset as usize - offset;

        if begin_offset >= end_offset {
            return cstr;
        }

        let operations =
            unsafe { std::slice::from_raw_parts(cigar.operations, cigar.max_operations as usize) };
        let mut last_op = operations[begin_offset];
        let mut last_op_length = 1;

        for i in 1..(end_offset - begin_offset) {
            let cur_op = operations[begin_offset + i];
            if cur_op == last_op {
                last_op_length += 1;
            } else {
                cstr.push_str(&format!("{}", last_op_length));
                cstr.push(last_op as u8 as char);
                last_op = cur_op;
                last_op_length = 1;
            }
        }
        cstr.push_str(&format!("{}", last_op_length));
        cstr.push(last_op as u8 as char);
        cstr
    }

    pub fn cigar_wfa(&self) -> Vec<u8> {
        if AlignmentScope::from(self.attributes.inner.alignment_scope) == AlignmentScope::Score {
            return Vec::new();
        }

        let cigar = unsafe { (*self.inner).cigar.as_ref() }
            .expect("CIGAR is null, alignment might have failed or scope was Score");

        if cigar.operations.is_null() || cigar.begin_offset > cigar.end_offset {
            return Vec::new();
        }

        let cigar_str = unsafe {
            let begin_offset = cigar.begin_offset;
            let cigar_operations = cigar.operations.offset(begin_offset as isize) as *const u8;
            let cigar_length = (cigar.end_offset - begin_offset) as usize;
            slice::from_raw_parts(cigar_operations, cigar_length)
        };
        cigar_str.to_vec()
    }

    pub fn matching(
        &self,
        pattern: &[u8],
        text: &[u8],
        flank_len: Option<usize>,
    ) -> (String, String, String) {
        let offset = flank_len.unwrap_or(0);

        let mut pattern_iter = pattern.iter().peekable();
        let mut text_iter = text.iter().peekable();

        if offset > 0 {
            text_iter.nth(offset - 1);
            pattern_iter.nth(offset - 1);
        }

        let mut pattern_alg = String::new();
        let mut ops_alg = String::new();
        let mut text_alg = String::new();

        let cigar = unsafe { (*self.inner).cigar.as_ref() }.unwrap();
        let operations =
            unsafe { std::slice::from_raw_parts(cigar.operations, cigar.max_operations as usize) };

        let begin_offset = cigar.begin_offset as isize + offset as isize;
        let end_offset = cigar.end_offset as isize - offset as isize;

        for i in begin_offset..end_offset {
            match operations[i as usize] as u8 as char {
                'M' => {
                    if pattern_iter.peek() != text_iter.peek() {
                        ops_alg.push('X');
                    } else {
                        ops_alg.push('|');
                    }
                    pattern_alg.push(*pattern_iter.next().unwrap() as char);
                    text_alg.push(*text_iter.next().unwrap() as char);
                }
                'X' => {
                    if pattern_iter.peek() != text_iter.peek() {
                        ops_alg.push(' ');
                    } else {
                        ops_alg.push('X');
                    }
                    pattern_alg.push(*pattern_iter.next().unwrap() as char);
                    text_alg.push(*text_iter.next().unwrap() as char);
                }
                'I' => {
                    pattern_alg.push('-');
                    ops_alg.push(' ');
                    text_alg.push(*text_iter.next().unwrap() as char);
                }
                'D' => {
                    pattern_alg.push(*pattern_iter.next().unwrap() as char);
                    ops_alg.push(' ');
                    text_alg.push('-');
                }
                _ => panic!("Unknown cigar operation"),
            }
        }
        (pattern_alg, ops_alg, text_alg)
    }

    pub fn find_alignment_span(&self) -> ((usize, usize), (usize, usize)) {
        let mut pattern_start: usize = 0;
        let mut pattern_end: usize = 0;
        let mut text_start: usize = 0;
        let mut text_end: usize = 0;

        let mut pattern_index: usize = 0;
        let mut text_index: usize = 0;
        let mut is_span_started: bool = false;

        let cigar = unsafe { (*self.inner).cigar.as_ref() }.unwrap();
        let operations =
            unsafe { std::slice::from_raw_parts(cigar.operations, cigar.max_operations as usize) };

        for i in cigar.begin_offset..cigar.end_offset {
            match operations[i as usize] as u8 as char {
                'I' => text_index += 1,
                'D' => pattern_index += 1,
                'M' => {
                    if !is_span_started {
                        pattern_start = pattern_index;
                        text_start = text_index;
                        is_span_started = true;
                    }
                    pattern_index += 1;
                    text_index += 1;
                    pattern_end = pattern_index;
                    text_end = text_index;
                }
                'X' => {
                    pattern_index += 1;
                    text_index += 1;
                }
                _ => panic!("Unexpected operation"),
            }
        }
        ((pattern_start, pattern_end), (text_start, text_end))
    }

    pub fn set_heuristic(&mut self, heuristic: Heuristic) {
        unsafe {
            match heuristic {
                Heuristic::None => wfa2::wavefront_aligner_set_heuristic_none(self.inner),
                Heuristic::WFadaptive(
                    min_wavefront_length,
                    max_distance_threshold,
                    steps_between_cutoffs,
                ) => wfa2::wavefront_aligner_set_heuristic_wfadaptive(
                    self.inner,
                    min_wavefront_length,
                    max_distance_threshold,
                    steps_between_cutoffs,
                ),
                Heuristic::WFmash(
                    min_wavefront_length,
                    max_distance_threshold,
                    steps_between_cutoffs,
                ) => wfa2::wavefront_aligner_set_heuristic_wfmash(
                    self.inner,
                    min_wavefront_length,
                    max_distance_threshold,
                    steps_between_cutoffs,
                ),
                Heuristic::XDrop(xdrop, steps_between_cutoffs) => {
                    wfa2::wavefront_aligner_set_heuristic_xdrop(
                        self.inner,
                        xdrop,
                        steps_between_cutoffs,
                    )
                }
                Heuristic::ZDrop(zdrop, steps_between_cutoffs) => {
                    wfa2::wavefront_aligner_set_heuristic_zdrop(
                        self.inner,
                        zdrop,
                        steps_between_cutoffs,
                    )
                }
                Heuristic::BandedStatic(band_min_k, band_max_k) => {
                    wfa2::wavefront_aligner_set_heuristic_banded_static(
                        self.inner, band_min_k, band_max_k,
                    )
                }
                Heuristic::BandedAdaptive(band_min_k, band_max_k, steps_between_cutoffs) => {
                    wfa2::wavefront_aligner_set_heuristic_banded_adaptive(
                        self.inner,
                        band_min_k,
                        band_max_k,
                        steps_between_cutoffs,
                    )
                }
            }
        }
    }
}

pub fn cigar_wfa_to_ops(wfa_cigar: &[u8]) -> Vec<AlignmentOperation> {
    let mut operations = Vec::with_capacity(wfa_cigar.len()); // Pre-allocate

    for &byte_op in wfa_cigar {
        let op = match byte_op {
            b'M' => AlignmentOperation::Match,
            b'X' => AlignmentOperation::Subst,
            b'I' => AlignmentOperation::Ins,
            b'D' => AlignmentOperation::Del,
            _ => {
                panic!(
                    "Unexpected byte '{}' found in WFA CIGAR string.",
                    byte_op as char
                );
            }
        };
        operations.push(op);
    }

    operations
}

/// Indel Aligner (a.k.a Longest Common Subsequence - LCS)
pub struct WFAlignerIndel;

impl WFAlignerIndel {
    pub fn create_aligner(alignment_scope: AlignmentScope, memory_model: MemoryModel) -> WFAligner {
        let mut aligner = WFAligner::new(alignment_scope, memory_model);
        aligner.attributes = aligner.attributes.indel_penalties();
        unsafe {
            aligner.inner = wfa2::wavefront_aligner_new(&mut aligner.attributes.inner);
        }
        aligner
    }
}

/// Edit Aligner (a.k.a Levenshtein)
pub struct WFAlignerEdit;

impl WFAlignerEdit {
    pub fn create_aligner(alignment_scope: AlignmentScope, memory_model: MemoryModel) -> WFAligner {
        let mut aligner = WFAligner::new(alignment_scope, memory_model);
        aligner.attributes = aligner.attributes.edit_penalties();
        unsafe {
            aligner.inner = wfa2::wavefront_aligner_new(&mut aligner.attributes.inner);
        }
        aligner
    }
}

/// Gap-Linear Aligner (a.k.a Needleman-Wunsch)
pub struct WFAlignerGapLinear;

impl WFAlignerGapLinear {
    pub fn create_aligner_with_match(
        match_: i32,
        mismatch: i32,
        indel: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        let mut aligner = WFAligner::new(alignment_scope, memory_model);
        aligner.attributes = WFAttributes::default()
            .linear_penalties(match_, mismatch, indel)
            .alignment_scope(alignment_scope)
            .memory_model(memory_model);
        unsafe {
            aligner.inner = wfa2::wavefront_aligner_new(&mut aligner.attributes.inner);
        }
        aligner
    }

    pub fn create_aligner(
        mismatch: i32,
        indel: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        Self::create_aligner_with_match(0, mismatch, indel, alignment_scope, memory_model)
    }
}

/// Gap-Affine Aligner (a.k.a Smith-Waterman-Gotoh)
pub struct WFAlignerGapAffine;

impl WFAlignerGapAffine {
    pub fn create_aligner_with_match(
        match_: i32,
        mismatch: i32,
        gap_opening: i32,
        gap_extension: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        let mut aligner = WFAligner::new(alignment_scope, memory_model);
        aligner.attributes = WFAttributes::default()
            .affine_penalties(match_, mismatch, gap_opening, gap_extension)
            .alignment_scope(alignment_scope)
            .memory_model(memory_model);
        unsafe {
            aligner.inner = wfa2::wavefront_aligner_new(&mut aligner.attributes.inner);
        }
        aligner
    }

    pub fn create_aligner(
        mismatch: i32,
        gap_opening: i32,
        gap_extension: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        Self::create_aligner_with_match(
            0,
            mismatch,
            gap_opening,
            gap_extension,
            alignment_scope,
            memory_model,
        )
    }
}

/// Gap-Affine Dual-Cost Aligner (a.k.a. concave 2-pieces)
pub struct WFAlignerGapAffine2Pieces;

impl WFAlignerGapAffine2Pieces {
    #[allow(clippy::too_many_arguments)]
    pub fn create_aligner_with_match(
        match_: i32,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        let mut aligner = WFAligner::new(alignment_scope, memory_model);
        aligner.attributes = WFAttributes::default()
            .affine2p_penalties(
                match_,
                mismatch,
                gap_opening1,
                gap_extension1,
                gap_opening2,
                gap_extension2,
            )
            .alignment_scope(alignment_scope)
            .memory_model(memory_model);
        unsafe {
            aligner.inner = wfa2::wavefront_aligner_new(&mut aligner.attributes.inner);
        }
        aligner
    }

    pub fn create_aligner(
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        Self::create_aligner_with_match(
            0,
            mismatch,
            gap_opening1,
            gap_extension1,
            gap_opening2,
            gap_extension2,
            alignment_scope,
            memory_model,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PATTERN: &[u8] = b"AGCTAGTGTCAATGGCTACTTTTCAGGTCCT";
    const TEXT: &[u8] = b"AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";

    #[test]
    fn test_aligner_indel() {
        let mut aligner =
            WFAlignerIndel::create_aligner(AlignmentScope::Alignment, MemoryModel::MemoryHigh);
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 10);
        assert_eq!(aligner.cigar(None), "1M1I1D3M1I5M2I2D8M1I1M1I1M1I9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "A-GCTA-GTGTC--AATGGCTACT-T-T-TCAGGTCCT\n|  ||| |||||    |||||||| | | |||||||||\nAA-CTAAGTGTCGG--TGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_edit() {
        let mut aligner =
            WFAlignerEdit::create_aligner(AlignmentScope::Alignment, MemoryModel::MemoryHigh);
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 7);
        assert_eq!(aligner.cigar(None), "1M1X3M1I5M2X8M1I1M1I1M1I9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AGCTA-GTGTCAATGGCTACT-T-T-TCAGGTCCT\n| ||| |||||  |||||||| | | |||||||||\nAACTAAGTGTCGGTGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_gap_linear() {
        let mut aligner = WFAlignerGapLinear::create_aligner(
            6,
            2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -20);
        assert_eq!(aligner.cigar(None), "1M1I1D3M1I5M2I2D8M1I1M1I1M1I9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "A-GCTA-GTGTC--AATGGCTACT-T-T-TCAGGTCCT\n|  ||| |||||    |||||||| | | |||||||||\nAA-CTAAGTGTCGG--TGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_gap_affine() {
        let mut aligner = WFAlignerGapAffine::create_aligner(
            6,
            4,
            2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -40);
        assert_eq!(aligner.cigar(None), "1M1X3M1I5M2X8M3I1M1X9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AGCTA-GTGTCAATGGCTACT---TTTCAGGTCCT\n| ||| |||||  ||||||||   | |||||||||\nAACTAAGTGTCGGTGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_score_only() {
        let mut aligner = WFAlignerGapAffine::create_aligner(
            6,
            4,
            2,
            AlignmentScope::Score,
            MemoryModel::MemoryLow,
        );
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -40);
        assert_eq!(aligner.cigar(None), "");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(format!("{}\n{}\n{}", a, b, c), "\n\n");
    }

    #[test]
    fn test_aligner_gap_affine_2pieces() {
        let mut aligner = WFAlignerGapAffine2Pieces::create_aligner(
            6,
            2,
            2,
            4,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -34);
        assert_eq!(aligner.cigar(None), "1M1X3M1I5M2X8M1I1M1I1M1I9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AGCTA-GTGTCAATGGCTACT-T-T-TCAGGTCCT\n| ||| |||||  |||||||| | | |||||||||\nAACTAAGTGTCGGTGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_span_1() {
        let pattern = b"AATTTAAGTCTAGGCTACTTTC";
        let text = b"CCGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT";
        let mut aligner = WFAlignerGapAffine2Pieces::create_aligner(
            8,
            4,
            2,
            24,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        let status = aligner.align_ends_free(pattern, 0, 0, text, 0, text.len() as i32);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        let ((xstart, xend), (ystart, yend)) = aligner.find_alignment_span();
        assert_eq!(ystart, 13);
        assert_eq!(yend, 35);
        assert_eq!(xstart, 0);
        assert_eq!(xend, 22);
    }

    #[test]
    fn test_aligner_span_2() {
        let pattern = b"GGGATCCCCGAAAAAGCGGGTTTGGCAAAAGCAAATTTCCCGAGTAAGCAGGCAGAGATCGCGCCAGACGCTCCCCAGAGCAGGGCGTCATGCACAAGAAAGCTTTGCACTTTGCGAACCAACGATAGGTGGGGGTGCGTGGAGGATGGAACACGGACGGCCCGGCTTGCTGCCTTCCCAGGCCTGCAGTTTGCCCATCCACGTCAGGGCCTCAGCCTGGCCGAAAGAAAGAAATGGTCTGTGATCCCCC";
        let text = b"AGCAGGGCGTCATGCACAAGAAAGCTTTGCACTTTGCGAACCAACGATAGGTGGGGGTGCGTGGAGGATGGAACACGGACGGCCCGGCTTGCTGCCTTCCCAGGCCTGCAGTTTGCCCATCCACGTCAGGGCCTCAGCCTGGCCGAAAGAAAGAAATGGTCTGTGATCCCCCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATTCCCGGCTACAAGGACCCTTCGAGCCCCGTTCGCCGGCCGCGGACCCGGCCCCTCCCTCCCCGGCCGCTAGGGGGCGGGCCCGGATCACAGGACTGGAGCTGGGCGGAGACCCACGCTCGGAGCGGTTGTGAACTGGCAGGCGGTGGGCGCGGCTTCTGTGCCGTGCCCCGGGCACTCAGTCTTCCAACGGGGCCCCGGAGTCGAAGACAGTTCTAGGGTTCAGGGAGCGCGGGCGGCTCCTGGGCGGCGCCAGACTGCGGTGAGTTGGCCGGCGTGGGCCACCAACCCAATGCAGCCCAGGGCGGCGGCACGAGACAGAACAACGGCGAACAGGAGCAGGGAAAGCGCCTCCGATAGGCCAGGCCTAGGGACCTGCGGGGAGAGGGCGAGGTCAACACCCGGCATGGGCCTCTGATTGGCTCCTGGGACTCGCCCCGCCTACGCCCATAGGTGGGCCCGCACTCTTCCCTGCGCCCCGCCCCCGCCCCAACAGCCT";
        let mut aligner = WFAlignerGapAffine2Pieces::create_aligner(
            8,
            4,
            2,
            24,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        aligner.set_heuristic(Heuristic::None);
        let status = aligner.align_ends_free(pattern, 0, 0, text, 0, text.len() as i32);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        let ((xstart, xend), (ystart, yend)) = aligner.find_alignment_span();

        assert_eq!(ystart, 0);
        assert_eq!(yend, 172);
        assert_eq!(xstart, 78);
        assert_eq!(xend, 250);
    }

    #[test]
    fn test_aligner_ends_free_global() {
        let pattern = b"AATTTAAGTCTAGGCTACTTTC";
        let text = b"CCGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT";
        let mut aligner = WFAlignerGapAffine::create_aligner(
            6,
            4,
            2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        let status = aligner.align_ends_free(pattern, 0, 0, text, 0, text.len() as i32);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -36);
        assert_eq!(aligner.cigar(None), "13I9M1X12M15I");
        let (a, b, c) = aligner.matching(pattern, text, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "-------------AATTTAAGTCTAGGCTACTTTC---------------\n             ||||||||| ||||||||||||               \nCCGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT"
        );
    }

    #[test]
    fn test_aligner_ends_free_right_extent() {
        let pattern = b"AATTTAAGTCTGCTACTTTCACGCAGCT";
        let text = b"AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT";
        let mut aligner = WFAlignerGapAffine::create_aligner(
            6,
            4,
            2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        let status =
            aligner.align_ends_free(pattern, 0, pattern.len() as i32, text, 0, text.len() as i32);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -24);
        assert_eq!(aligner.cigar(None), "5M1X6M1I11M4D1M15I");
        let (a, b, c) = aligner.matching(pattern, text, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AATTTAAGTCTG-CTACTTTCACGCAGCT---------------\n||||| |||||| |||||||||||    |               \nAATTTCAGTCTGGCTACTTTCACG----TACGATGACAGACTCT"
        );
    }

    #[test]
    fn test_aligner_ends_free_left_extent() {
        let pattern = b"CTTTCACGTACGTGACAGTCTCT";
        let text = b"AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT";
        let mut aligner = WFAlignerGapAffine::create_aligner(
            6,
            4,
            2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        let status = aligner.align_ends_free(pattern, 0, 0, text, 0, 0);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -48);
        assert_eq!(aligner.cigar(None), "16I12M1I6M1X4M");
        let (a, b, c) = aligner.matching(pattern, text, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "----------------CTTTCACGTACG-TGACAGTCTCT\n                |||||||||||| |||||| ||||\nAATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT"
        );
    }

    #[test]
    fn test_aligner_ends_free_right_overlap() {
        let pattern = b"CGCGTCTGACTGACTGACTAAACTTTCATGTACCTGACA";
        let text = b"AAACTTTCACGTACGTGACATATAGCGATCGATGACT";
        let mut aligner = WFAlignerGapAffine::create_aligner(
            6,
            4,
            2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        let status = aligner.align_ends_free(pattern, 0, 0, text, 0, 0);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -92);
        assert_eq!(aligner.cigar(None), "19D9M1X4M1X5M17I");
        let (a, b, c) = aligner.matching(pattern, text, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "CGCGTCTGACTGACTGACTAAACTTTCATGTACCTGACA-----------------\n                   ||||||||| |||| |||||                 \n-------------------AAACTTTCACGTACGTGACATATAGCGATCGATGACT"
        );
    }

    #[test]
    fn test_clipping_score() {
        let text_lf = b"AAGGAGCTGAGAATTGTTCTTCCAGATACCTTTCCGACCTCTTCTTGGTT";
        let text_rf = b"GGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCCGCATCCTGGGTT";

        let pattern_lf = b"AAGGAGCTGAGAATTGTTCGTCCAGATACCTTTCCGACCTCTTCTTGGTT";
        let pattern_rf = b"GGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCTGCATCCTGGGTT";

        let motif = b"ATTT";

        let text = [text_lf, &motif.repeat(10)[..], text_rf].concat();
        let pattern = [pattern_lf, &motif.repeat(8)[..], pattern_rf].concat();

        let mut aligner = WFAlignerGapAffine2Pieces::create_aligner(
            8,
            4,
            2,
            24,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        let status = aligner.align_end_to_end(&pattern, &text);

        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -36);
        assert_eq!(aligner.cigar(None), "19M1X62M8I37M1X12M");
        let (a, b, c) = aligner.matching(&pattern, &text, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AAGGAGCTGAGAATTGTTCGTCCAGATACCTTTCCGACCTCTTCTTGGTTATTTATTTATTTATTTATTTATTTATTTATTT--------GGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCTGCATCCTGGGTT\n||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||        ||||||||||||||||||||||||||||||||||||| ||||||||||||\nAAGGAGCTGAGAATTGTTCTTCCAGATACCTTTCCGACCTCTTCTTGGTTATTTATTTATTTATTTATTTATTTATTTATTTATTTATTTGGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCCGCATCCTGGGTT"
        );
        assert_eq!(aligner.cigar_score(), -36);
        assert_eq!(aligner.cigar_score_clipped(50), -20);
        assert_eq!(aligner.cigar(Some(50)), "32M8I");
        let (a, b, c) = aligner.matching(&pattern, &text, Some(50));
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "ATTTATTTATTTATTTATTTATTTATTTATTT--------\n||||||||||||||||||||||||||||||||        \nATTTATTTATTTATTTATTTATTTATTTATTTATTTATTT"
        );
    }

    #[test]
    fn test_memory_modes() {
        let expected_cigar = "1M1X3M1I5M2X8M3I1M1X9M";
        let expected_matching = "AGCTA-GTGTCAATGGCTACT---TTTCAGGTCCT\n| ||| |||||  ||||||||   | |||||||||\nAACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";
        let expected_score = -48;

        struct Test {
            memory_mode: MemoryModel,
        }

        let tests = vec![
            Test {
                memory_mode: MemoryModel::MemoryHigh,
            },
            Test {
                memory_mode: MemoryModel::MemoryMed,
            },
            Test {
                memory_mode: MemoryModel::MemoryLow,
            },
            // Test {
            //     memory_mode: MemoryModel::MemoryUltraLow,
            // },
        ];

        for test in tests {
            let mut aligner = WFAlignerGapAffine2Pieces::create_aligner(
                8,
                4,
                2,
                24,
                1,
                AlignmentScope::Alignment,
                test.memory_mode,
            );
            let status = aligner.align_end_to_end(PATTERN, TEXT);
            assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
            assert_eq!(aligner.score(), expected_score);
            assert_eq!(aligner.cigar_score(), expected_score);
            assert_eq!(aligner.cigar_score_clipped(0), expected_score);
            assert_eq!(aligner.cigar(None), expected_cigar);
            let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
            assert_eq!(format!("{}\n{}\n{}", a, b, c), expected_matching);
        }
    }

    #[test]
    fn test_set_heuristic() {
        let mut aligner = WFAlignerGapAffine::create_aligner(
            6,
            4,
            2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        aligner.set_heuristic(Heuristic::WFmash(1, 2, 3));
        aligner.set_heuristic(Heuristic::BandedStatic(1, 2));
        aligner.set_heuristic(Heuristic::BandedAdaptive(1, 2, 3));
        aligner.set_heuristic(Heuristic::WFadaptive(1, 2, 3));
        aligner.set_heuristic(Heuristic::XDrop(1, 2));
        aligner.set_heuristic(Heuristic::ZDrop(1, 2));
        aligner.set_heuristic(Heuristic::None);
    }

    #[test]
    fn test_invalid_sequence() {
        let read = b"GCTGCTACTGGGGTGTCCCCTCTCAAAGGACAAACCCAGGATCTACAGATGTGTGTGCTAAGCCATGTATGCACATGCACGTGTGTGTGTATATATTTAACCTATCTGTATATATGTATTATGTAAACATGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCCTGCTGGCATATCTGACTATAACTGACCACCTCACAGTCCATTCTGATCTCTATATATGTATTATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATTATGTAAACATGAGTTCCCTGCTGGCATATCTGATTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATTATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCTATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGATTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATTATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGATCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGGCTGGCATATCTGATTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGATTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCCCGCTGGCTTTTCCATGACTTCCTTATCCAGCTGTGAGAACCCTGACTCTTACTACCCATACTGTATTGACTTATTT";
        let allele = b"GCTGCTACTGGGGTGTCCCCTCTCAAAGGACAAACCCAGGATCTACAGATGTGTGTGCTAAGCCATGTATGCACACGCACGTGTGTGTGTATATATTTAACCTATCTGTATATATGTATTATGTAAACATGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGACTTCCTACTGGCATATCTGACTGTAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGATTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTTCATTCCGATCTGTATATAAGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTGTAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATGTATGTATCATGTAAACACGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCCGATCTGTATATAAGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTGTAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGCATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGTCTATAACCGACCACCTTAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGTCCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGTCTATAACCGACCACCTTAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGCATATATGTATAATATATATTATATATGGTCCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGTCTATAACCGACCACCTTAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCCCGCTGGCTTTTCCATGACTTCCTTATCCAGCTGTGAGAACCCTGACTCTTACTACTGTATTGACTTATTTGTGAAACCT";

        let mut aligner = WFAlignerGapAffine2Pieces::create_aligner(
            8,
            4,
            2,
            24,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryUltraLow,
        );
        let _status = aligner.align_end_to_end(read, allele);
        assert_eq!(_status, AlignmentStatus::StatusUnattainable);
        assert_eq!(aligner.score(), -2147483648);

        aligner.set_heuristic(Heuristic::None);
        let _status = aligner.align_end_to_end(read, allele);
        assert_eq!(_status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -881);
    }

    #[test]
    fn test_get_penalties() {
        let aligner =
            WFAlignerEdit::create_aligner(AlignmentScope::Alignment, MemoryModel::MemoryLow);
        assert_eq!(aligner.get_penalties(), (0, 4, Some(2), None, None, None));

        let aligner =
            WFAlignerIndel::create_aligner(AlignmentScope::Alignment, MemoryModel::MemoryLow);
        assert_eq!(aligner.get_penalties(), (0, 4, Some(2), None, None, None));

        let aligner = WFAlignerGapLinear::create_aligner(
            12,
            24,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        assert_eq!(aligner.get_penalties(), (0, 12, Some(24), None, None, None));

        let aligner = WFAlignerGapAffine::create_aligner(
            12,
            24,
            2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        assert_eq!(
            aligner.get_penalties(),
            (0, 12, Some(24), Some(2), None, None)
        );

        let aligner = WFAlignerGapAffine2Pieces::create_aligner(
            12,
            24,
            2,
            48,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        assert_eq!(
            aligner.get_penalties(),
            (0, 12, Some(24), Some(2), Some(48), Some(1))
        );
    }

    #[test]
    fn test_get_and_decode_sam_cigar() {
        let pattern = b"TCTTTACTCTT";
        let text = b"TCTTTACTCTT";
        let mut aligner = WFAlignerGapAffine::create_aligner(
            4,
            6,
            2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );

        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);

        let sam_cigar_buffer = aligner.get_sam_cigar(true);
        assert!(
            !sam_cigar_buffer.is_empty(),
            "SAM CIGAR buffer should not be empty"
        );

        let decoded_cigar = WFAligner::decode_sam_cigar(&sam_cigar_buffer);

        // Expected result for identical sequences (11 matches)
        // The raw buffer encodes length << 4 | op_code. '=' is op_code 7.
        // So, 11= should be encoded as (11 << 4) | 7 = 176 | 7 = 183
        let expected_raw_buffer = vec![183]; // 11=
        assert_eq!(
            sam_cigar_buffer, expected_raw_buffer,
            "Raw SAM CIGAR buffer mismatch"
        );

        let expected_decoded_cigar = vec![(11, '=')]; // 11 matches ('=' because show_mismatches=true)
        assert_eq!(
            decoded_cigar, expected_decoded_cigar,
            "Decoded SAM CIGAR mismatch"
        );

        // Test with show_mismatches = false (should use 'M')
        let sam_cigar_buffer_m = aligner.get_sam_cigar(false);
        // 'M' is op_code 0. (11 << 4) | 0 = 176
        let expected_raw_buffer_m = vec![176]; // 11M
        assert_eq!(
            sam_cigar_buffer_m, expected_raw_buffer_m,
            "Raw SAM CIGAR buffer mismatch (M)"
        );

        let decoded_cigar_m = WFAligner::decode_sam_cigar(&sam_cigar_buffer_m);
        let expected_decoded_cigar_m = vec![(11, 'M')]; // 11 matches ('M')
        assert_eq!(
            decoded_cigar_m, expected_decoded_cigar_m,
            "Decoded SAM CIGAR mismatch (M)"
        );

        // Test with a simple mismatch
        let pattern_diff = b"TCTTTACTCTT";
        let text_diff = b"TCTTTACTATT";
        let mut aligner_diff = WFAlignerGapAffine::create_aligner(
            4,
            6,
            2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );
        let status_diff = aligner_diff.align_end_to_end(pattern_diff, text_diff);
        assert_eq!(status_diff, AlignmentStatus::StatusAlgCompleted);

        let sam_cigar_buffer_diff = aligner_diff.get_sam_cigar(true);

        let expected_raw_diff = vec![135, 24, 39];
        assert_eq!(
            sam_cigar_buffer_diff, expected_raw_diff,
            "Raw SAM CIGAR buffer mismatch (diff)"
        );

        let decoded_cigar_diff = WFAligner::decode_sam_cigar(&sam_cigar_buffer_diff);
        let expected_decoded_diff = vec![(8, '='), (1, 'X'), (2, '=')];
        assert_eq!(
            decoded_cigar_diff, expected_decoded_diff,
            "Decoded SAM CIGAR mismatch (diff)"
        );

        // Test with show_mismatches = false
        let sam_cigar_buffer_diff_m = aligner_diff.get_sam_cigar(false);
        // Expected: 11M => (11<<4)|0 = 176
        let expected_raw_diff_m = vec![176];
        assert_eq!(
            sam_cigar_buffer_diff_m, expected_raw_diff_m,
            "Raw SAM CIGAR buffer mismatch (diff, M)"
        );

        let decoded_cigar_diff_m = WFAligner::decode_sam_cigar(&sam_cigar_buffer_diff_m);
        let expected_decoded_diff_m = vec![(11, 'M')];
        assert_eq!(
            decoded_cigar_diff_m, expected_decoded_diff_m,
            "Decoded SAM CIGAR mismatch (diff, M)"
        );
    }
}
