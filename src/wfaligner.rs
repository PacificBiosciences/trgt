use crate::wfa2;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryModel {
    MemoryHigh,
    MemoryMed,
    MemoryLow,
    MemoryUltraLow,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WfaOp {
    Match,
    Subst,
    Ins,
    Del,
}

impl WfaOp {
    fn from_u8(op_char: u8) -> Self {
        match op_char {
            b'M' => WfaOp::Match,
            b'X' => WfaOp::Subst,
            b'I' => WfaOp::Ins,
            b'D' => WfaOp::Del,
            _ => panic!("Invalid alignment operation character {}", op_char),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WfaAlign {
    /// WFA alignment score (cost).
    pub score: i32,
    /// Start position of alignment in the reference (text). 0-based.
    pub ystart: usize,
    /// Start position of alignment in the query (pattern). 0-based.
    pub xstart: usize,
    /// End position of alignment in the reference (text). 0-based, exclusive.
    pub yend: usize,
    /// End position of alignment in the query (pattern). 0-based, exclusive.
    pub xend: usize,
    /// Length of the reference sequence (text) involved in the alignment.
    pub ylen: usize,
    /// Length of the query sequence (pattern) involved in the alignment.
    pub xlen: usize,
    /// Vector of alignment operations.
    pub operations: Vec<WfaOp>,
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

#[derive(Debug, PartialEq, Eq)]
pub enum Penalties {
    Indel, // Conceptually: mismatch=inf, indel=1
    Edit,  // Conceptually: mismatch=1, indel=1
    Linear {
        match_: i32,
        mismatch: i32,
        indel: i32,
    },
    Affine {
        match_: i32,
        mismatch: i32,
        gap_opening: i32,
        gap_extension: i32,
    },
    Affine2p {
        match_: i32,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
    },
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

impl fmt::Display for AlignmentStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlignmentStatus::StatusAlgCompleted => write!(f, "StatusAlgCompleted"),
            AlignmentStatus::StatusAlgPartial => write!(f, "StatusAlgPartial"),
            AlignmentStatus::StatusMaxStepsReached => write!(f, "StatusMaxStepsReached"),
            AlignmentStatus::StatusOOM => write!(f, "StatusOOM"),
            AlignmentStatus::StatusUnattainable => write!(f, "StatusUnattainable"),
        }
    }
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

pub struct WFAlignerBuilder {
    attributes: WFAttributes,
    penalty_set: bool,
    heuristic: Option<Heuristic>,
}

impl WFAlignerBuilder {
    pub fn new(alignment_scope: AlignmentScope, memory_model: MemoryModel) -> Self {
        let attributes = WFAttributes::default()
            .memory_model(memory_model)
            .alignment_scope(alignment_scope);
        Self {
            attributes,
            penalty_set: false,
            heuristic: None,
        }
    }

    /// Configure for indel penalties (Longest Common Subsequence - LCS)
    pub fn indel(mut self) -> Self {
        self.attributes = self.attributes.indel_penalties();
        self.penalty_set = true;
        self
    }

    /// Configure for edit penalties (Levenshtein)
    pub fn edit(mut self) -> Self {
        self.attributes = self.attributes.edit_penalties();
        self.penalty_set = true;
        self
    }

    /// Configure for gap-linear penalties (Needleman-Wunsch) with match_ = 0
    pub fn linear(self, mismatch: i32, indel: i32) -> Self {
        self.linear_with_match(0, mismatch, indel)
    }

    /// Configure for gap-linear penalties (Needleman-Wunsch) with explicit match score
    pub fn linear_with_match(mut self, match_: i32, mismatch: i32, indel: i32) -> Self {
        self.attributes = self.attributes.linear_penalties(match_, mismatch, indel);
        self.penalty_set = true;
        self
    }

    /// Configure for gap-affine penalties (Smith-Waterman-Gotoh) with match_ = 0
    pub fn affine(self, mismatch: i32, gap_opening: i32, gap_extension: i32) -> Self {
        self.affine_with_match(0, mismatch, gap_opening, gap_extension)
    }

    /// Configure for gap-affine penalties (Smith-Waterman-Gotoh) with explicit match score
    pub fn affine_with_match(
        mut self,
        match_: i32,
        mismatch: i32,
        gap_opening: i32,
        gap_extension: i32,
    ) -> Self {
        self.attributes =
            self.attributes
                .affine_penalties(match_, mismatch, gap_opening, gap_extension);
        self.penalty_set = true;
        self
    }

    /// Configure for gap-affine dual-cost penalties (concave 2-pieces) with match_ = 0
    pub fn affine2p(
        self,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
    ) -> Self {
        self.affine2p_with_match(
            0,
            mismatch,
            gap_opening1,
            gap_extension1,
            gap_opening2,
            gap_extension2,
        )
    }

    /// Configure for gap-affine dual-cost penalties (concave 2-pieces) with explicit match score
    #[allow(clippy::too_many_arguments)]
    pub fn affine2p_with_match(
        mut self,
        match_: i32,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
    ) -> Self {
        self.attributes = self.attributes.affine2p_penalties(
            match_,
            mismatch,
            gap_opening1,
            gap_extension1,
            gap_opening2,
            gap_extension2,
        );
        self.penalty_set = true;
        self
    }

    /// Set a heuristic for the aligner
    pub fn with_heuristic(mut self, heuristic: Heuristic) -> Self {
        self.heuristic = Some(heuristic);
        self
    }

    /// Build the WFAligner with the configured settings
    pub fn build(self) -> WFAligner {
        if !self.penalty_set {
            panic!("Must set a penalty model before building the aligner");
        }

        let mut aligner = WFAligner {
            attributes: self.attributes,
            inner: std::ptr::null_mut(),
        };

        unsafe {
            aligner.inner = wfa2::wavefront_aligner_new(&mut aligner.attributes.inner);
        }

        if let Some(heuristic) = self.heuristic {
            aligner.set_heuristic(heuristic);
        }

        aligner
    }
}

// TODO: Unify different Cigar wrappers
/// Represents a single operation: (length, op).
pub type CigarOp = (usize, char);

pub struct WFAligner {
    attributes: WFAttributes,
    inner: *mut wfa2::wavefront_aligner_t,
}

impl WFAligner {
    /// Create a builder for configuring a WFAligner
    pub fn builder(alignment_scope: AlignmentScope, memory_model: MemoryModel) -> WFAlignerBuilder {
        WFAlignerBuilder::new(alignment_scope, memory_model)
    }

    pub fn get_penalties(&self) -> Penalties {
        match DistanceMetric::from(self.attributes.inner.distance_metric) {
            DistanceMetric::Indel => Penalties::Indel,
            DistanceMetric::Edit => Penalties::Edit,
            DistanceMetric::GapLinear => Penalties::Linear {
                match_: self.attributes.inner.linear_penalties.match_,
                mismatch: self.attributes.inner.linear_penalties.mismatch,
                indel: self.attributes.inner.linear_penalties.indel,
            },
            DistanceMetric::GapAffine => Penalties::Affine {
                match_: self.attributes.inner.affine_penalties.match_,
                mismatch: self.attributes.inner.affine_penalties.mismatch,
                gap_opening: self.attributes.inner.affine_penalties.gap_opening,
                gap_extension: self.attributes.inner.affine_penalties.gap_extension,
            },
            DistanceMetric::GapAffine2p => Penalties::Affine2p {
                match_: self.attributes.inner.affine2p_penalties.match_,
                mismatch: self.attributes.inner.affine2p_penalties.mismatch,
                gap_opening1: self.attributes.inner.affine2p_penalties.gap_opening1,
                gap_extension1: self.attributes.inner.affine2p_penalties.gap_extension1,
                gap_opening2: self.attributes.inner.affine2p_penalties.gap_opening2,
                gap_extension2: self.attributes.inner.affine2p_penalties.gap_extension2,
            },
        }
    }

    pub fn get_heuristics(&self) -> Heuristic {
        let h = &self.attributes.inner.heuristic;
        match h.strategy {
            wfa2::wf_heuristic_strategy_wf_heuristic_none => Heuristic::None,
            wfa2::wf_heuristic_strategy_wf_heuristic_banded_static => {
                Heuristic::BandedStatic(h.min_k, h.max_k)
            }
            wfa2::wf_heuristic_strategy_wf_heuristic_banded_adaptive => {
                Heuristic::BandedAdaptive(h.min_k, h.max_k, h.steps_between_cutoffs)
            }
            wfa2::wf_heuristic_strategy_wf_heuristic_wfadaptive => Heuristic::WFadaptive(
                h.min_wavefront_length,
                h.max_distance_threshold,
                h.steps_between_cutoffs,
            ),
            wfa2::wf_heuristic_strategy_wf_heuristic_wfmash => Heuristic::WFmash(
                h.min_wavefront_length,
                h.max_distance_threshold,
                h.steps_between_cutoffs,
            ),
            wfa2::wf_heuristic_strategy_wf_heuristic_xdrop => {
                Heuristic::XDrop(h.xdrop, h.steps_between_cutoffs)
            }
            wfa2::wf_heuristic_strategy_wf_heuristic_zdrop => {
                Heuristic::ZDrop(h.zdrop, h.steps_between_cutoffs)
            }
            _ => panic!("Unknown heuristic strategy: {}", h.strategy),
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
    fn set_alignment_end_to_end(&mut self) {
        unsafe {
            wfa2::wavefront_aligner_set_alignment_end_to_end(self.inner);
        }
    }

    fn set_alignment_ends_free(
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

    fn cigar_score_edit_get_operations_score(operation: char, op_length: i32) -> i32 {
        match operation {
            'M' => 0,                     // Match has zero cost in edit distance
            'X' | 'D' | 'I' => op_length, // Mismatch, Deletion, Insertion cost 1 each
            _ => panic!("Invalid operation: {}", operation),
        }
    }

    // Note: Indel distance (LCS) is usually about maximizing matches, not minimizing penalties.
    // For scoring purposes, we treat it like Edit distance in a cost minimization framework
    fn cigar_score_indel_get_operations_score(operation: char, op_length: i32) -> i32 {
        match operation {
            'M' => 0,                     // Match has zero cost
            'X' | 'D' | 'I' => op_length, // Non-matches cost 1 each
            _ => panic!("Invalid operation: {}", operation),
        }
    }

    fn cigar_score_gap_linear_get_operations_score(
        operation: char,
        op_length: i32,
        penalties: &wfa2::linear_penalties_t,
    ) -> i32 {
        match operation {
            'M' => op_length * penalties.match_,
            'X' => op_length * penalties.mismatch,
            'D' | 'I' => op_length * penalties.indel,
            _ => panic!("Invalid operation: {}", operation),
        }
    }

    fn cigar_score_gap_affine_get_operations_score(
        operation: char,
        op_length: i32,
        penalties: &wfa2::affine_penalties_t,
    ) -> i32 {
        match operation {
            'M' => op_length * penalties.match_,
            'X' => op_length * penalties.mismatch,
            'D' | 'I' => penalties.gap_opening + penalties.gap_extension * op_length,
            _ => panic!("Invalid operation: {}", operation),
        }
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

    pub fn cigar_score_clipped(&self, flank_len: usize) -> i32 {
        if AlignmentScope::from(self.attributes.inner.alignment_scope) == AlignmentScope::Score {
            panic!("Cannot clip when AlignmentScope is Score");
        }

        let cigar = unsafe { (*self.inner).cigar.as_ref() }.unwrap();
        let begin_offset = cigar.begin_offset as isize + flank_len as isize;
        // Ensure end_offset doesn't go below begin_offset
        let end_offset =
            std::cmp::max(begin_offset, cigar.end_offset as isize - flank_len as isize);

        if begin_offset >= end_offset {
            return 0; // No operations in the clipped region
        }

        let operations: Vec<char> =
            unsafe { std::slice::from_raw_parts(cigar.operations, cigar.max_operations as usize) }
                .iter()
                .map(|&op| op as u8 as char)
                .collect();

        let mut score = 0;
        let mut op_length = 0;
        // Initialize last_op with the first operation in the clipped range
        let mut last_op: Option<char> = Some(operations[begin_offset as usize]);

        for i in begin_offset..end_offset {
            let cur_op = operations[i as usize];
            if Some(cur_op) != last_op {
                // Calculate score for the completed segment of the previous operation
                if let Some(op) = last_op {
                    match DistanceMetric::from(self.attributes.inner.distance_metric) {
                        DistanceMetric::Indel => {
                            // Add cost for Indel (positive value)
                            score += Self::cigar_score_indel_get_operations_score(op, op_length);
                        }
                        DistanceMetric::Edit => {
                            // Add cost for Edit (positive value)
                            score += Self::cigar_score_edit_get_operations_score(op, op_length);
                        }
                        DistanceMetric::GapLinear => {
                            // Subtract penalty for gap-linear (negative value)
                            score -= Self::cigar_score_gap_linear_get_operations_score(
                                op,
                                op_length,
                                &self.attributes.inner.linear_penalties,
                            );
                        }
                        DistanceMetric::GapAffine => {
                            // Subtract penalty for gap-affine (negative value)
                            score -= Self::cigar_score_gap_affine_get_operations_score(
                                op,
                                op_length,
                                &self.attributes.inner.affine_penalties,
                            );
                        }
                        DistanceMetric::GapAffine2p => {
                            // Subtract penalty for gap-affine-2p (negative value)
                            score -= Self::cigar_score_gap_affine2p_get_operations_score(
                                op,
                                op_length,
                                &self.attributes.inner.affine2p_penalties,
                            );
                        }
                    };
                }
                op_length = 0; // Reset length for the new operation
            }
            last_op = Some(cur_op);
            op_length += 1;
        }

        // Add score for the last segment
        if let Some(op) = last_op {
            match DistanceMetric::from(self.attributes.inner.distance_metric) {
                DistanceMetric::Indel => {
                    // Add cost for Indel (positive value)
                    score += Self::cigar_score_indel_get_operations_score(op, op_length);
                }
                DistanceMetric::Edit => {
                    // Add cost for Edit (positive value)
                    score += Self::cigar_score_edit_get_operations_score(op, op_length);
                }
                DistanceMetric::GapLinear => {
                    // Subtract penalty for gap-linear (negative value)
                    score -= Self::cigar_score_gap_linear_get_operations_score(
                        op,
                        op_length,
                        &self.attributes.inner.linear_penalties,
                    );
                }
                DistanceMetric::GapAffine => {
                    // Subtract penalty for gap-affine (negative value)
                    score -= Self::cigar_score_gap_affine_get_operations_score(
                        op,
                        op_length,
                        &self.attributes.inner.affine_penalties,
                    );
                }
                DistanceMetric::GapAffine2p => {
                    // Subtract penalty for gap-affine-2p (negative value)
                    score -= Self::cigar_score_gap_affine2p_get_operations_score(
                        op,
                        op_length,
                        &self.attributes.inner.affine2p_penalties,
                    );
                }
            };
        }
        score
    }

    pub fn set_heuristic(&mut self, heuristic: Heuristic) {
        match heuristic {
            Heuristic::None => {
                self.attributes.inner.heuristic.strategy =
                    wfa2::wf_heuristic_strategy_wf_heuristic_none;
                unsafe {
                    wfa2::wavefront_aligner_set_heuristic_none(self.inner);
                }
            }
            Heuristic::WFadaptive(min_len, max_dist, steps) => {
                self.attributes.inner.heuristic.strategy =
                    wfa2::wf_heuristic_strategy_wf_heuristic_wfadaptive;
                self.attributes.inner.heuristic.min_wavefront_length = min_len;
                self.attributes.inner.heuristic.max_distance_threshold = max_dist;
                self.attributes.inner.heuristic.steps_between_cutoffs = steps;
                unsafe {
                    wfa2::wavefront_aligner_set_heuristic_wfadaptive(
                        self.inner, min_len, max_dist, steps,
                    );
                }
            }
            Heuristic::WFmash(min_len, max_dist, steps) => {
                self.attributes.inner.heuristic.strategy =
                    wfa2::wf_heuristic_strategy_wf_heuristic_wfmash;
                self.attributes.inner.heuristic.min_wavefront_length = min_len;
                self.attributes.inner.heuristic.max_distance_threshold = max_dist;
                self.attributes.inner.heuristic.steps_between_cutoffs = steps;
                unsafe {
                    wfa2::wavefront_aligner_set_heuristic_wfmash(
                        self.inner, min_len, max_dist, steps,
                    );
                }
            }
            Heuristic::XDrop(xdrop, steps) => {
                self.attributes.inner.heuristic.strategy =
                    wfa2::wf_heuristic_strategy_wf_heuristic_xdrop;
                self.attributes.inner.heuristic.xdrop = xdrop;
                self.attributes.inner.heuristic.steps_between_cutoffs = steps;
                unsafe {
                    wfa2::wavefront_aligner_set_heuristic_xdrop(self.inner, xdrop, steps);
                }
            }
            Heuristic::ZDrop(zdrop, steps) => {
                self.attributes.inner.heuristic.strategy =
                    wfa2::wf_heuristic_strategy_wf_heuristic_zdrop;
                self.attributes.inner.heuristic.zdrop = zdrop;
                self.attributes.inner.heuristic.steps_between_cutoffs = steps;
                unsafe {
                    wfa2::wavefront_aligner_set_heuristic_zdrop(self.inner, zdrop, steps);
                }
            }
            Heuristic::BandedStatic(min_k, max_k) => {
                self.attributes.inner.heuristic.strategy =
                    wfa2::wf_heuristic_strategy_wf_heuristic_banded_static;
                self.attributes.inner.heuristic.min_k = min_k;
                self.attributes.inner.heuristic.max_k = max_k;
                unsafe {
                    wfa2::wavefront_aligner_set_heuristic_banded_static(self.inner, min_k, max_k);
                }
            }
            Heuristic::BandedAdaptive(min_k, max_k, steps) => {
                self.attributes.inner.heuristic.strategy =
                    wfa2::wf_heuristic_strategy_wf_heuristic_banded_adaptive;
                self.attributes.inner.heuristic.min_k = min_k;
                self.attributes.inner.heuristic.max_k = max_k;
                self.attributes.inner.heuristic.steps_between_cutoffs = steps;
                unsafe {
                    wfa2::wavefront_aligner_set_heuristic_banded_adaptive(
                        self.inner, min_k, max_k, steps,
                    );
                }
            }
        }
    }

    // TODO: Refactor
    // TODO: THIS DOES NOT WORK WITH BIWFA
    pub fn get_alignment(&self) -> WfaAlign {
        let cigar = unsafe { (*self.inner).cigar.as_ref() }.unwrap();
        let raw_operations =
            unsafe { std::slice::from_raw_parts(cigar.operations, cigar.max_operations as usize) };

        let begin = cigar.begin_offset as usize;
        let end = cigar.end_offset as usize;

        let mut operations = Vec::with_capacity(raw_operations.len());
        for op in &raw_operations[begin..end] {
            let operation = WfaOp::from_u8(*op as u8);
            operations.push(operation);
        }

        let pattern_len = unsafe { (*self.inner).sequences.pattern_length } as usize;
        let text_len = unsafe { (*self.inner).sequences.text_length } as usize;

        let (xstart, xend, ystart, yend);

        // Check if alignment was end-to-end
        let is_global = unsafe {
            (*self.inner).alignment_form.span == wfa2::alignment_span_t_alignment_end2end
        };
        // For global alignment, span is always the full sequence lengths
        if is_global {
            xstart = 0;
            ystart = 0;
            xend = pattern_len;
            yend = text_len;
        } else {
            // For ends-free, we need to calculate the span
            let mut pattern_index = 0;
            let mut text_index = 0;
            let mut is_span_started = false;
            let (mut current_pattern_start, mut current_text_start) = (0, 0);
            let (mut current_pattern_end, mut current_text_end) = (0, 0);

            for op in &operations {
                match op {
                    WfaOp::Ins => {
                        text_index += 1;
                    }
                    WfaOp::Del => {
                        pattern_index += 1;
                    }
                    WfaOp::Match | WfaOp::Subst => {
                        if !is_span_started {
                            current_pattern_start = pattern_index;
                            current_text_start = text_index;
                            is_span_started = true;
                        }
                        pattern_index += 1;
                        text_index += 1;
                        current_pattern_end = pattern_index;
                        current_text_end = text_index;
                    }
                }
            }

            // Handle cases where alignment might be purely indels (unlikely for ends-free) or if the loop didn't run. If span never started, ends remain 0. If it did start, the last match/mismatch determined the end
            xstart = current_pattern_start;
            ystart = current_text_start;
            xend = current_pattern_end;
            yend = current_text_end;
        }

        WfaAlign {
            score: cigar.score,
            ystart,
            yend,
            xstart,
            xend,
            ylen: text_len,
            xlen: pattern_len,
            operations,
        }
    }

    // TODO: Deprecate and use get_alignment?
    // TODO: THIS DOES NOT WORK WITH BIWFA
    pub fn get_alignment_span(&self) -> ((usize, usize), (usize, usize)) {
        // Check if alignment was end-to-end
        let is_global = unsafe {
            (*self.inner).alignment_form.span == wfa2::alignment_span_t_alignment_end2end
        };

        if is_global {
            let pattern_len = unsafe { (*self.inner).sequences.pattern_length } as usize;
            let text_len = unsafe { (*self.inner).sequences.text_length } as usize;
            ((0, pattern_len), (0, text_len))
        } else {
            let cigar = unsafe { (*self.inner).cigar.as_ref() }.unwrap();
            let raw_operations = unsafe {
                std::slice::from_raw_parts(cigar.operations, cigar.max_operations as usize)
            };

            let begin = cigar.begin_offset as usize;
            let end = cigar.end_offset as usize;

            let mut pattern_index = 0;
            let mut text_index = 0;
            let mut is_span_started = false;
            let (mut pattern_start, mut pattern_end, mut text_start, mut text_end) = (0, 0, 0, 0);

            for &op in &raw_operations[begin..end] {
                match op as u8 {
                    b'I' => text_index += 1,
                    b'D' => pattern_index += 1,
                    b'M' | b'X' => {
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
                    _ => panic!("Unexpected operation"),
                }
            }
            ((pattern_start, pattern_end), (text_start, text_end))
        }
    }

    pub fn cigar_operations(&self) -> Vec<u8> {
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
            std::slice::from_raw_parts(cigar_operations, cigar_length)
        };
        cigar_str.to_vec()
    }

    // TODO: Possibly remove safety checks?
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
                    std::slice::from_raw_parts(sam_cigar_buffer_ptr, sam_cigar_length as usize);
                cigar_buffer_slice.to_vec()
            } else {
                Vec::new()
            }
        }
    }

    pub fn decode_sam_cigar(sam_cigar_buffer: &[u32]) -> Vec<CigarOp> {
        const SAM_CIGAR_LEN_SHIFT: u32 = 4;
        const SAM_CIGAR_OP_MASK: u32 = 0xF;
        sam_cigar_buffer
            .iter()
            .map(|&encoded_op| {
                let len = encoded_op >> SAM_CIGAR_LEN_SHIFT; // Length is in the upper 28 bits
                let op_code = encoded_op & SAM_CIGAR_OP_MASK; // Operation code is in the lower 4 bits
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
                (len as usize, op_char)
            })
            .collect()
    }

    /// Counts the number of match ('M') operations in the CIGAR string.
    // TODO: Possibly remove safety checks?
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

    #[allow(dead_code)]
    fn cigar_string(&self, flank_len: Option<usize>) -> String {
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

    #[allow(dead_code)]
    fn matching(
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use WfaOp::*;

    const PATTERN: &[u8] = b"AGCTAGTGTCAATGGCTACTTTTCAGGTCCT";
    const TEXT: &[u8] = b"AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";

    #[test]
    fn test_aligner_indel() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .indel()
            .build();
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 10);
        assert_eq!(aligner.cigar_string(None), "1M1I1D3M1I5M2I2D8M1I1M1I1M1I9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "A-GCTA-GTGTC--AATGGCTACT-T-T-TCAGGTCCT\n|  ||| |||||    |||||||| | | |||||||||\nAA-CTAAGTGTCGG--TGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_edit() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .edit()
            .build();
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 7);
        assert_eq!(aligner.cigar_string(None), "1M1X3M1I5M2X8M1I1M1I1M1I9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AGCTA-GTGTCAATGGCTACT-T-T-TCAGGTCCT\n| ||| |||||  |||||||| | | |||||||||\nAACTAAGTGTCGGTGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_gap_linear() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .linear(6, 2)
            .build();
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -20);
        assert_eq!(aligner.cigar_string(None), "1M1I1D3M1I5M2I2D8M1I1M1I1M1I9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "A-GCTA-GTGTC--AATGGCTACT-T-T-TCAGGTCCT\n|  ||| |||||    |||||||| | | |||||||||\nAA-CTAAGTGTCGG--TGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_gap_affine() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine(6, 4, 2)
            .build();
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -40);
        assert_eq!(aligner.cigar_string(None), "1M1X3M1I5M2X8M3I1M1X9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AGCTA-GTGTCAATGGCTACT---TTTCAGGTCCT\n| ||| |||||  ||||||||   | |||||||||\nAACTAAGTGTCGGTGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_score_only() {
        let mut aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryLow)
            .affine(6, 4, 2)
            .build();
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -40);
        assert_eq!(aligner.cigar_string(None), "");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(format!("{}\n{}\n{}", a, b, c), "\n\n");
    }

    #[test]
    fn test_aligner_gap_affine_2pieces() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine2p(6, 2, 2, 4, 1)
            .build();
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -34);
        assert_eq!(aligner.cigar_string(None), "1M1X3M1I5M2X8M1I1M1I1M1I9M");
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
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine2p(8, 4, 2, 24, 1)
            .build();
        let status = aligner.align_ends_free(pattern, 0, 0, text, 0, text.len() as i32);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        let ((xstart, xend), (ystart, yend)) = aligner.get_alignment_span();
        assert_eq!(ystart, 13);
        assert_eq!(yend, 35);
        assert_eq!(xstart, 0);
        assert_eq!(xend, 22);
    }

    #[test]
    fn test_aligner_span_2() {
        let pattern = b"GGGATCCCCGAAAAAGCGGGTTTGGCAAAAGCAAATTTCCCGAGTAAGCAGGCAGAGATCGCGCCAGACGCTCCCCAGAGCAGGGCGTCATGCACAAGAAAGCTTTGCACTTTGCGAACCAACGATAGGTGGGGGTGCGTGGAGGATGGAACACGGACGGCCCGGCTTGCTGCCTTCCCAGGCCTGCAGTTTGCCCATCCACGTCAGGGCCTCAGCCTGGCCGAAAGAAAGAAATGGTCTGTGATCCCCC";
        let text = b"AGCAGGGCGTCATGCACAAGAAAGCTTTGCACTTTGCGAACCAACGATAGGTGGGGGTGCGTGGAGGATGGAACACGGACGGCCCGGCTTGCTGCCTTCCCAGGCCTGCAGTTTGCCCATCCACGTCAGGGCCTCAGCCTGGCCGAAAGAAAGAAATGGTCTGTGATCCCCCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATTCCCGGCTACAAGGACCCTTCGAGCCCCGTTCGCCGGCCGCGGACCCGGCCCCTCCCTCCCCGGCCGCTAGGGGGCGGGCCCGGATCACAGGACTGGAGCTGGGCGGAGACCCACGCTCGGAGCGGTTGTGAACTGGCAGGCGGTGGGCGCGGCTTCTGTGCCGTGCCCCGGGCACTCAGTCTTCCAACGGGGCCCCGGAGTCGAAGACAGTTCTAGGGTTCAGGGAGCGCGGGCGGCTCCTGGGCGGCGCCAGACTGCGGTGAGTTGGCCGGCGTGGGCCACCAACCCAATGCAGCCCAGGGCGGCGGCACGAGACAGAACAACGGCGAACAGGAGCAGGGAAAGCGCCTCCGATAGGCCAGGCCTAGGGACCTGCGGGGAGAGGGCGAGGTCAACACCCGGCATGGGCCTCTGATTGGCTCCTGGGACTCGCCCCGCCTACGCCCATAGGTGGGCCCGCACTCTTCCCTGCGCCCCGCCCCCGCCCCAACAGCCT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine2p(8, 4, 2, 24, 1)
            .build();
        aligner.set_heuristic(Heuristic::None);
        let status = aligner.align_ends_free(pattern, 0, 0, text, 0, text.len() as i32);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        let ((xstart, xend), (ystart, yend)) = aligner.get_alignment_span();

        assert_eq!(ystart, 0);
        assert_eq!(yend, 172);
        assert_eq!(xstart, 78);
        assert_eq!(xend, 250);
    }

    #[test]
    fn test_aligner_ends_free_global() {
        let pattern = b"AATTTAAGTCTAGGCTACTTTC";
        let text = b"CCGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();
        let status = aligner.align_ends_free(pattern, 0, 0, text, 0, text.len() as i32);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -36);
        assert_eq!(aligner.cigar_string(None), "13I9M1X12M15I");
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
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();
        let status =
            aligner.align_ends_free(pattern, 0, pattern.len() as i32, text, 0, text.len() as i32);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -24);
        assert_eq!(aligner.cigar_string(None), "5M1X6M1I11M4D1M15I");
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
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();
        let status = aligner.align_ends_free(pattern, 0, 0, text, 0, 0);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -48);
        assert_eq!(aligner.cigar_string(None), "16I12M1I6M1X4M");
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
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();
        let status = aligner.align_ends_free(pattern, 0, 0, text, 0, 0);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -92);
        assert_eq!(aligner.cigar_string(None), "19D9M1X4M1X5M17I");
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

        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine2p(8, 4, 2, 24, 1)
            .build();
        let status = aligner.align_end_to_end(&pattern, &text);

        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -36);
        assert_eq!(aligner.cigar_string(None), "19M1X62M8I37M1X12M");
        let (a, b, c) = aligner.matching(&pattern, &text, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AAGGAGCTGAGAATTGTTCGTCCAGATACCTTTCCGACCTCTTCTTGGTTATTTATTTATTTATTTATTTATTTATTTATTT--------GGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCTGCATCCTGGGTT\n||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||        ||||||||||||||||||||||||||||||||||||| ||||||||||||\nAAGGAGCTGAGAATTGTTCTTCCAGATACCTTTCCGACCTCTTCTTGGTTATTTATTTATTTATTTATTTATTTATTTATTTATTTATTTGGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCCGCATCCTGGGTT"
        );
        assert_eq!(aligner.cigar_score(), -36);
        assert_eq!(aligner.cigar_score_clipped(50), -20);
        assert_eq!(aligner.cigar_string(Some(50)), "32M8I");
        let (a, b, c) = aligner.matching(&pattern, &text, Some(50));
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "ATTTATTTATTTATTTATTTATTTATTTATTT--------\n||||||||||||||||||||||||||||||||        \nATTTATTTATTTATTTATTTATTTATTTATTTATTTATTT"
        );

        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .indel()
            .with_heuristic(Heuristic::None)
            .build();
        let status = aligner.align_end_to_end(&pattern, &text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 12);
        assert_eq!(aligner.cigar_score(), 12);
        assert_eq!(aligner.cigar_score_clipped(19), 10);
        assert_eq!(aligner.cigar_score_clipped(0), 12);
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
            let mut aligner = WFAligner::builder(AlignmentScope::Alignment, test.memory_mode)
                .affine2p(8, 4, 2, 24, 1)
                .build();
            let status = aligner.align_end_to_end(PATTERN, TEXT);
            assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
            assert_eq!(aligner.score(), expected_score);
            assert_eq!(aligner.cigar_score(), expected_score);
            assert_eq!(aligner.cigar_score_clipped(0), expected_score);
            assert_eq!(aligner.cigar_string(None), expected_cigar);
            let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
            assert_eq!(format!("{}\n{}\n{}", a, b, c), expected_matching);
        }
    }

    #[test]
    fn test_set_heuristic() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();
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

        let mut aligner =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
                .affine2p(8, 4, 2, 24, 1)
                .build();
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
        let aligner_edit = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .edit()
            .build();
        assert_eq!(aligner_edit.get_penalties(), Penalties::Edit);

        let aligner_indel = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .indel()
            .build();
        assert_eq!(aligner_indel.get_penalties(), Penalties::Indel);

        let aligner_linear = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .linear(12, 24)
            .build();
        assert_eq!(
            aligner_linear.get_penalties(),
            Penalties::Linear {
                match_: 0,
                mismatch: 12,
                indel: 24
            }
        );

        let aligner_affine = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine(12, 24, 2)
            .build();
        assert_eq!(
            aligner_affine.get_penalties(),
            Penalties::Affine {
                match_: 0,
                mismatch: 12,
                gap_opening: 24,
                gap_extension: 2
            }
        );

        let aligner_affine2p =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .affine2p(12, 24, 2, 48, 1)
                .build();
        assert_eq!(
            aligner_affine2p.get_penalties(),
            Penalties::Affine2p {
                match_: 0,
                mismatch: 12,
                gap_opening1: 24,
                gap_extension1: 2,
                gap_opening2: 48,
                gap_extension2: 1
            }
        );

        let aligner_affine_match =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .affine_with_match(-5, 12, 24, 2)
                .build();
        assert_eq!(
            aligner_affine_match.get_penalties(),
            Penalties::Affine {
                match_: -5,
                mismatch: 12,
                gap_opening: 24,
                gap_extension: 2
            }
        );
    }

    #[test]
    fn test_builder_pattern() {
        let aligner_edit = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .edit()
            .build();
        assert_eq!(aligner_edit.get_penalties(), Penalties::Edit);

        let aligner_affine = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine(12, 24, 2)
            .with_heuristic(Heuristic::WFadaptive(10, 50, 100))
            .build();
        assert_eq!(
            aligner_affine.get_penalties(),
            Penalties::Affine {
                match_: 0,
                mismatch: 12,
                gap_opening: 24,
                gap_extension: 2
            }
        );

        let aligner_affine_heuristic =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .affine(12, 24, 2)
                .with_heuristic(Heuristic::WFadaptive(10, 50, 100))
                .build();
        assert_eq!(
            aligner_affine_heuristic.get_penalties(),
            Penalties::Affine {
                match_: 0,
                mismatch: 12,
                gap_opening: 24,
                gap_extension: 2
            }
        );

        let aligner_affine2p =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .affine2p(12, 24, 2, 48, 1)
                .build();
        assert_eq!(
            aligner_affine2p.get_penalties(),
            Penalties::Affine2p {
                match_: 0,
                mismatch: 12,
                gap_opening1: 24,
                gap_extension1: 2,
                gap_opening2: 48,
                gap_extension2: 1
            }
        );

        let aligner_linear = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .linear_with_match(-5, 12, 24)
            .build();
        assert_eq!(
            aligner_linear.get_penalties(),
            Penalties::Linear {
                match_: -5,
                mismatch: 12,
                indel: 24
            }
        );
    }

    #[test]
    fn test_get_and_decode_sam_cigar() {
        let pattern = b"TCTTTACTCTT";
        let text = b"TCTTTACTCTT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(4, 6, 2)
            .build();

        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);

        let sam_cigar_buffer = aligner.get_sam_cigar(true);
        assert!(
            !sam_cigar_buffer.is_empty(),
            "SAM CIGAR buffer should not be empty"
        );

        let decoded_cigar = WFAligner::decode_sam_cigar(&sam_cigar_buffer);

        // Expected result for identical sequences (11 matches), The raw buffer encodes length << 4 | op_code. '=' is op_code 7. So, 11= should be encoded as (11 << 4) | 7 = 176 | 7 = 183
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

        // Test with show_mismatches = false
        let sam_cigar_buffer_m = aligner.get_sam_cigar(false);
        // 'M' is op_code 0. (11 << 4) | 0 = 176
        let expected_raw_buffer_m = vec![176]; // 11M
        assert_eq!(
            sam_cigar_buffer_m, expected_raw_buffer_m,
            "Raw SAM CIGAR buffer mismatch (M)"
        );

        let decoded_cigar_m = WFAligner::decode_sam_cigar(&sam_cigar_buffer_m);
        let expected_decoded_cigar_m: Vec<CigarOp> = vec![(11, 'M')]; // 11 matches ('M')
        assert_eq!(
            decoded_cigar_m, expected_decoded_cigar_m,
            "Decoded SAM CIGAR mismatch (M)"
        );

        let pattern_diff = b"TCTTTACTCTT";
        let text_diff = b"TCTTTACTATT";
        let mut aligner_diff =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .affine(4, 6, 2)
                .build();
        let status_diff = aligner_diff.align_end_to_end(pattern_diff, text_diff);
        assert_eq!(status_diff, AlignmentStatus::StatusAlgCompleted);

        let sam_cigar_buffer_diff = aligner_diff.get_sam_cigar(true);

        let expected_raw_diff = vec![135, 24, 39];
        assert_eq!(
            sam_cigar_buffer_diff, expected_raw_diff,
            "Raw SAM CIGAR buffer mismatch (diff)"
        );

        let decoded_cigar_diff = WFAligner::decode_sam_cigar(&sam_cigar_buffer_diff);
        let expected_decoded_diff: Vec<CigarOp> = vec![(8, '='), (1, 'X'), (2, '=')];
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
        let expected_decoded_diff_m: Vec<CigarOp> = vec![(11, 'M')];
        assert_eq!(
            decoded_cigar_diff_m, expected_decoded_diff_m,
            "Decoded SAM CIGAR mismatch (diff, M)"
        );
    }

    #[test]
    fn test_get_heuristics() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(12, 24, 2)
            .with_heuristic(Heuristic::None)
            .build();
        assert_eq!(aligner.get_heuristics(), Heuristic::None);

        let wf_adaptive_heuristic = Heuristic::WFadaptive(5, 25, 50);
        aligner.set_heuristic(wf_adaptive_heuristic);
        assert_eq!(aligner.get_heuristics(), wf_adaptive_heuristic);

        let banded_static_heuristic = Heuristic::BandedStatic(5, 20);
        aligner.set_heuristic(banded_static_heuristic);
        assert_eq!(aligner.get_heuristics(), banded_static_heuristic);

        let xdrop_heuristic = Heuristic::XDrop(15, 5);
        aligner.set_heuristic(xdrop_heuristic);
        assert_eq!(aligner.get_heuristics(), xdrop_heuristic);

        let aligner_with_heuristic =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
                .affine(12, 24, 2)
                .with_heuristic(Heuristic::WFadaptive(10, 50, 100))
                .build();
        assert_eq!(
            aligner_with_heuristic.get_heuristics(),
            Heuristic::WFadaptive(10, 50, 100)
        );
        assert_eq!(
            aligner_with_heuristic.get_penalties(),
            Penalties::Affine {
                match_: 0,
                mismatch: 12,
                gap_opening: 24,
                gap_extension: 2
            }
        );
    }

    #[test]
    fn test_get_alignment_global() {
        let pattern = b"AGCTAGTGTCAATGGCTACTTTTCAGGTCCT";
        let text = b"AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(1, 5, 1)
            .build();
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);

        let alignment = aligner.get_alignment();
        assert_eq!(aligner.score(), -18);

        let expected_ops = vec![
            Match, Subst, Match, Match, Match, Ins, Match, Match, Match, Match, Match, Subst,
            Subst, Match, Match, Match, Match, Match, Match, Match, Match, Ins, Ins, Ins, Match,
            Subst, Match, Match, Match, Match, Match, Match, Match, Match, Match,
        ];

        assert_eq!(alignment.score, -18);
        assert_eq!(alignment.xlen, pattern.len());
        assert_eq!(alignment.ylen, text.len());
        assert_eq!(alignment.operations, expected_ops);

        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, 31);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 35);

        let ((xstart, xend), (ystart, yend)) = aligner.get_alignment_span();
        assert_eq!(alignment.xstart, xstart);
        assert_eq!(alignment.xend, xend);
        assert_eq!(alignment.ystart, ystart);
        assert_eq!(alignment.yend, yend);
    }

    // #[test]
    // fn test_get_alignment_biwfa_global() {
    //     let pattern = b"AGCTAGTGTCAATGGCTACTTTTCAGGTCCT";
    //     let text = b"AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";
    //     let mut aligner =
    //         WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
    //             .affine(1, 5, 1)
    //             .build();
    //     let status = aligner.align_end_to_end(pattern, text);
    //     assert_eq!(status, AlignmentStatus::StatusAlgCompleted);

    //     let alignment = aligner.get_alignment();
    //     assert_eq!(aligner.score(), -2147483648);
    //     assert_eq!(aligner.cigar_score(), -18);

    //     let expected_ops = vec![
    //         Match, Subst, Match, Match, Match, Ins, Match, Match, Match, Match, Match, Subst,
    //         Subst, Match, Match, Match, Match, Match, Match, Match, Match, Ins, Ins, Ins, Match,
    //         Subst, Match, Match, Match, Match, Match, Match, Match, Match, Match,
    //     ];

    //     // WARN: score, xlen, and ylen cannot be relied on in biWFA mode!
    //     assert_eq!(alignment.score, -2147483648);
    //     assert_ne!(alignment.xlen, pattern.len());
    //     assert_ne!(alignment.ylen, text.len());

    //     assert_eq!(alignment.operations, expected_ops);

    //     assert_eq!(alignment.xstart, 0);
    //     assert_eq!(alignment.xend, 31);
    //     assert_eq!(alignment.ystart, 0);
    //     assert_eq!(alignment.yend, 35);

    //     let ((xstart, xend), (ystart, yend)) = aligner.get_alignment_span();
    //     assert_eq!(alignment.xstart, xstart);
    //     assert_eq!(alignment.xend, xend);
    //     assert_eq!(alignment.ystart, ystart);
    //     assert_eq!(alignment.yend, yend);
    // }

    #[test]
    fn test_get_alignment_ends_free() {
        let pattern = b"AGTGTCAATGGCTAC";
        let text = b"GGGGGGGGGGAGTGTCAATGGCTACGGGGGGGGGG";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(1, 5, 1)
            .build();
        let status =
            aligner.align_ends_free(pattern, 0, 0, text, text.len() as i32, text.len() as i32);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        let alignment = aligner.get_alignment();
        assert_eq!(aligner.score(), 0);

        let expected_ops = vec![
            Ins, Ins, Ins, Ins, Ins, Ins, Ins, Ins, Ins, Ins, Match, Match, Match, Match, Match,
            Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Ins, Ins, Ins,
            Ins, Ins, Ins, Ins, Ins, Ins, Ins,
        ];

        assert_eq!(alignment.score, 0);
        assert_eq!(alignment.xlen, pattern.len());
        assert_eq!(alignment.ylen, text.len());
        assert_eq!(alignment.operations, expected_ops);

        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, pattern.len());
        assert_eq!(alignment.ystart, 10);
        assert_eq!(alignment.yend, 25);

        let ((xstart, xend), (ystart, yend)) = aligner.get_alignment_span();
        assert_eq!(alignment.xstart, xstart);
        assert_eq!(alignment.xend, xend);
        assert_eq!(alignment.ystart, ystart);
        assert_eq!(alignment.yend, yend);
    }
}
