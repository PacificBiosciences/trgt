use crate::utils::Result;

pub fn create_writer<T, F>(output_prefix: &str, output_suffix: &str, f: F) -> Result<T>
where
    F: FnOnce(&str) -> Result<T>,
{
    let output_path = format!("{}.{}", output_prefix, output_suffix);
    f(&output_path)
}
