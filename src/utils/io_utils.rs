use crate::utils::Result;
use std::path::Path;

pub fn create_writer<T, F>(output_prefix: &Path, output_suffix: &str, f: F) -> Result<T>
where
    F: FnOnce(&str) -> Result<T>,
{
    let mut output_path = output_prefix.to_path_buf().into_os_string();
    output_path.push(format!(".{output_suffix}"));
    f(output_path.to_str().unwrap())
}
