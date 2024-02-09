mod builder;
mod events;
mod hmm;
mod purity;
mod spans;
mod utils;

pub use builder::build_hmm;
pub use events::get_events;
pub use events::HmmEvent;
pub use hmm::Hmm;
pub use purity::calc_purity;
pub use spans::*;
pub use utils::collapse_labels;
pub use utils::count_motifs;
pub use utils::replace_invalid_bases;
