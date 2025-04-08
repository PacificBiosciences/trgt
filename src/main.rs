use clap::Parser;
use std::time;
use trgt::{
    cli::{init_verbose, Cli, Command, FULL_VERSION},
    commands::{genotype, merge, plot, validate},
    utils::{handle_error_and_exit, Result},
};

fn runner() -> Result<()> {
    let cli = Cli::parse();
    init_verbose(&cli);
    log::info!(
        "Running {}-{} [{}]",
        env!("CARGO_PKG_NAME"),
        *FULL_VERSION,
        cli.command.name()
    );

    let start_timer = time::Instant::now();
    match cli.command {
        Command::Genotype(args) => {
            log::trace!("Genotype arguments: {:#?}", args);
            genotype::trgt(args)?
        }
        Command::Plot(args) => {
            log::trace!("Plot arguments: {:#?}", args);
            plot::trvz(args)?
        }
        Command::Validate(args) => {
            log::trace!("Validate arguments: {:#?}", args);
            validate::validate(args)?
        }
        Command::Merge(args) => {
            log::trace!("Merge arguments: {:#?}", args);
            merge::merge(args)?
        }
    }
    log::info!("Total execution time: {:.2?}", start_timer.elapsed());
    log::info!("{} end", env!("CARGO_PKG_NAME"));
    Ok(())
}

fn main() {
    if let Err(e) = runner() {
        handle_error_and_exit(e);
    }
}
