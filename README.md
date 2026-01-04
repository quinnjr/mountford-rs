# mountford-rs

Rust implementation of the Mountford plugin for PluMA - Dissimilarity index computation based on presence/absence data.

## Overview

This plugin computes the Mountford dissimilarity index (Mountford, 1962) between ecological samples. The Mountford index is a presence/absence based measure that emphasizes species shared between samples.

### Mountford Index

The Mountford similarity index is defined as:

```
M_AB = 2J / (2ab - (a+b)J)
```

Where:
- `a` = number of species in sample A
- `b` = number of species in sample B
- `J` = number of species common to both samples

The plugin outputs dissimilarity (1 - similarity), where:
- 0 = identical species composition
- 1 = no shared species

## Installation

```bash
cargo install --path .
```

## Usage

### Command Line

```bash
mountford abundance.csv dissimilarity.csv
```

### Input Format

CSV with species as rows and samples as columns:

```csv
species,Sample1,Sample2,Sample3
Bacteroides,150,0,200
Firmicutes,100,250,0
Proteobacteria,0,100,150
```

Values > 0 indicate presence; values = 0 indicate absence.

### Output Format

CSV dissimilarity matrix:

```csv
,Sample1,Sample2,Sample3
Sample1,0.000000,0.800000,0.500000
Sample2,0.800000,0.000000,0.666667
Sample3,0.500000,0.666667,0.000000
```

### As a Library

```rust
use mountford_rs::MountfordPlugin;

let mut plugin = MountfordPlugin::new();
plugin.input("abundance.csv")?;
plugin.run();

let diss = plugin.dissimilarity();
println!("Dissimilarity[0][1]: {}", diss[0][1]);

plugin.output("output.csv")?;
```

## Features

- Presence/absence based dissimilarity
- Parallel computation for large sample sets
- CSV input/output compatible with R's vegan package
- Symmetric dissimilarity matrix

## Building

```bash
cargo build --release
```

## Testing

```bash
cargo test
```

## Benchmarking

```bash
cargo bench
```

## References

- Mountford, M.D. (1962). An index of similarity and its application to classificatory problems. In: P.W. Murphy (ed.), Progress in Soil Zoology, 43-50.
- Original PluMA plugin: https://github.com/movingpictures83/Mountford

## License

MIT
