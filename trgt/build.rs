use std::error::Error;
use vergen::EmitBuilder;

fn main() -> Result<(), Box<dyn Error>> {
    match EmitBuilder::builder()
        .fail_on_error()
        .custom_build_rs(".") // override such that it will re-run whenever we have a file change in this folder
        .all_git()
        .git_describe(true, false, Some("ThisPatternShouldNotMatchAnythingEver"))
        .emit()
    {
        Ok(_) => {}
        Err(_e) => {
            println!("cargo:rustc-env=VERGEN_GIT_DESCRIBE=unknown");
        }
    }
    Ok(())
}
