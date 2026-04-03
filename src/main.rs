use clap::Parser;
pub use montecarlo::{
    error::{Error, Result},
    run,
};

/// Monte carlo simulation of photons absorption by skin
#[derive(Parser, Debug)]
#[command(about, long_about = None)]
struct Args {
    /// Number of 1'000s of photon to simulate
    #[arg(short, long, default_value_t = 100)]
    kphoton: u16,

    /// Wavelength of the photons, in nm (many even ìnteger between 600 and 1000 allowed)
    #[arg(short, long, default_value_t = 700)]
    wavelenght: u16,

    /// Verbosity
    #[arg(short, long, default_value_t = false)]
    verbose: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    run(args.kphoton, args.wavelenght, args.verbose)
}
