use clap::Parser;
use trgt::{
    cli::{init_verbose, Cli, Command, FULL_VERSION},
    commands::{genotype, merge, plot, validate},
    utils::{handle_error_and_exit, Result},
};

fn runner() -> Result<()> {
    let cli = Cli::parse();
    init_verbose(&cli);
    let subcommand_name = match cli.command {
        Command::Genotype(_) => "genotype",
        Command::Plot(_) => "plot",
        Command::Validate(_) => "validate",
        Command::Merge(_) => "merge",
    };

    log::info!(
        "Running {}-{} [{}]",
        env!("CARGO_PKG_NAME"),
        *FULL_VERSION,
        subcommand_name
    );
    match cli.command {
        Command::Genotype(args) => genotype::trgt(args)?,
        Command::Plot(args) => plot::trvz(args)?,
        Command::Validate(args) => validate::validate(args)?,
        Command::Merge(args) => merge::merge(args)?,
    }
    log::info!("{} end", env!("CARGO_PKG_NAME"));
    Ok(())
}

fn main() {
    if let Err(e) = runner() {
        handle_error_and_exit(e);
    }
}
