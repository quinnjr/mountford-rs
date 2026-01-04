//! Benchmarks for Mountford using Criterion v0.8.0

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use mountford_rs::MountfordPlugin;
use std::hint::black_box;

/// Generate a random-ish abundance matrix
fn generate_matrix(
    num_samples: usize,
    num_species: usize,
) -> (Vec<String>, Vec<String>, Vec<Vec<f64>>) {
    let samples: Vec<String> = (0..num_samples)
        .map(|i| format!("Sample{}", i + 1))
        .collect();
    let species: Vec<String> = (0..num_species)
        .map(|i| format!("Species{}", i + 1))
        .collect();

    let abundance: Vec<Vec<f64>> = (0..num_samples)
        .map(|i| {
            (0..num_species)
                .map(|j| {
                    // Deterministic "random" pattern
                    if (i * 7 + j * 13) % 3 == 0 {
                        ((i + j) % 100) as f64
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();

    (samples, species, abundance)
}

fn bench_dissimilarity_computation(c: &mut Criterion) {
    let mut group = c.benchmark_group("dissimilarity_computation");

    for (num_samples, num_species) in [(10, 50), (25, 100), (50, 200), (100, 500)] {
        let (samples, species, abundance) = generate_matrix(num_samples, num_species);
        let plugin = MountfordPlugin::from_matrix(samples, species, abundance);

        group.bench_with_input(
            BenchmarkId::new("compute", format!("{}x{}", num_samples, num_species)),
            &plugin,
            |b, p| {
                b.iter(|| {
                    let mut p = p.clone();
                    p.run();
                    black_box(p.dissimilarity().len())
                })
            },
        );
    }

    group.finish();
}

fn bench_large_matrices(c: &mut Criterion) {
    let mut group = c.benchmark_group("large_matrices");
    group.sample_size(20);

    for (num_samples, num_species) in [(200, 1000), (500, 500)] {
        let (samples, species, abundance) = generate_matrix(num_samples, num_species);
        let plugin = MountfordPlugin::from_matrix(samples, species, abundance);

        group.bench_with_input(
            BenchmarkId::new("compute", format!("{}x{}", num_samples, num_species)),
            &plugin,
            |b, p| {
                b.iter(|| {
                    let mut p = p.clone();
                    p.run();
                    black_box(p.dissimilarity().len())
                })
            },
        );
    }

    group.finish();
}

fn bench_full_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("full_pipeline");

    for (num_samples, num_species) in [(25, 100), (50, 200)] {
        let (samples, species, abundance) = generate_matrix(num_samples, num_species);

        group.bench_with_input(
            BenchmarkId::new(
                "from_matrix_and_run",
                format!("{}x{}", num_samples, num_species),
            ),
            &(samples.clone(), species.clone(), abundance.clone()),
            |b, (s, sp, a)| {
                b.iter(|| {
                    let mut plugin = MountfordPlugin::from_matrix(s.clone(), sp.clone(), a.clone());
                    plugin.run();
                    black_box(plugin.dissimilarity().len())
                })
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_dissimilarity_computation,
    bench_large_matrices,
    bench_full_pipeline
);
criterion_main!(benches);
