mod builder;
mod events;
pub mod hmm_model;
mod purity;
pub mod spans;
pub mod utils;

pub use builder::build_hmm;
pub use events::{get_events, HmmEvent};
pub use hmm_model::Hmm;
pub use purity::calc_purity;
pub use spans::Annotation;
pub use utils::{collapse_labels, count_motifs, replace_invalid_bases};
