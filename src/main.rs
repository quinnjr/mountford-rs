//! Mountford CLI - Dissimilarity index computation tool
//!
//! Usage: mountford <input.csv> <output.csv>

use mountford_rs::MountfordPlugin;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 3 {
        eprintln!("Usage: {} <input.csv> <output.csv>", args[0]);
        eprintln!();
        eprintln!("Computes Mountford dissimilarity index between samples.");
        eprintln!("Input: CSV with species (rows) x samples (columns)");
        eprintln!("Output: CSV dissimilarity matrix");
        process::exit(1);
    }

    let input_file = &args[1];
    let output_file = &args[2];

    let mut plugin = MountfordPlugin::new();

    // Input phase
    if let Err(e) = plugin.input(input_file) {
        eprintln!("Error reading input file '{}': {}", input_file, e);
        process::exit(1);
    }

    eprintln!(
        "Loaded {} samples x {} species",
        plugin.num_samples(),
        plugin.num_species()
    );

    // Run phase
    plugin.run();

    let n = plugin.num_samples();
    let pairs = n * (n - 1) / 2;
    eprintln!("Computed {} pairwise dissimilarities", pairs);

    // Output phase
    if let Err(e) = plugin.output(output_file) {
        eprintln!("Error writing output file '{}': {}", output_file, e);
        process::exit(1);
    }

    eprintln!("Results written to '{}'", output_file);
}
