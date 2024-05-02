pub type Result<T> = std::result::Result<T, String>;

pub fn handle_error_and_exit(err: String) -> ! {
    log::error!("{}", err);
    std::process::exit(1);
}
