//! Mountford - Dissimilarity Index computation
//!
//! This plugin computes the Mountford dissimilarity index (Mountford, 1962)
//! between samples based on presence/absence of species.
//!
//! The Mountford index is defined as:
//! M_AB = 2J / (2ab - (a+b)J)
//!
//! Where:
//! - a = number of species in sample A
//! - b = number of species in sample B
//! - J = number of species common to both samples

use rayon::prelude::*;
use std::path::Path;

/// Mountford plugin for PluMA
#[derive(Clone)]
pub struct MountfordPlugin {
    /// Sample names
    samples: Vec<String>,
    /// Species names
    species: Vec<String>,
    /// Abundance matrix (samples x species)
    abundance: Vec<Vec<f64>>,
    /// Presence/absence matrix (samples x species)
    presence: Vec<Vec<bool>>,
    /// Dissimilarity matrix (samples x samples)
    dissimilarity: Vec<Vec<f64>>,
}

impl MountfordPlugin {
    /// Create a new empty plugin
    pub fn new() -> Self {
        Self {
            samples: Vec::new(),
            species: Vec::new(),
            abundance: Vec::new(),
            presence: Vec::new(),
            dissimilarity: Vec::new(),
        }
    }

    /// Load abundance matrix from CSV file
    /// Format: First column is species names, subsequent columns are sample abundances
    pub fn input<P: AsRef<Path>>(&mut self, path: P) -> Result<(), Box<dyn std::error::Error>> {
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .flexible(true)
            .from_path(path)?;

        // Get sample names from headers (skip first column which is species name)
        let headers = reader.headers()?.clone();
        self.samples = headers.iter().skip(1).map(|s| s.to_string()).collect();

        let num_samples = self.samples.len();

        // Read abundance data
        self.species.clear();
        self.abundance = vec![vec![]; num_samples];

        for result in reader.records() {
            let record = result?;
            let mut iter = record.iter();

            // First column is species name
            if let Some(species_name) = iter.next() {
                self.species.push(species_name.to_string());

                // Remaining columns are abundances for each sample
                for (i, field) in iter.enumerate() {
                    if i < num_samples {
                        let value: f64 = field.trim().parse().unwrap_or(0.0);
                        self.abundance[i].push(value);
                    }
                }
            }
        }

        // Convert abundance to presence/absence
        self.compute_presence();

        Ok(())
    }

    /// Convert abundance matrix to presence/absence
    fn compute_presence(&mut self) {
        self.presence = self
            .abundance
            .iter()
            .map(|sample| sample.iter().map(|&v| v > 0.0).collect())
            .collect();
    }

    /// Compute the Mountford dissimilarity matrix
    pub fn run(&mut self) {
        if self.presence.is_empty() {
            return;
        }

        let n = self.samples.len();
        self.dissimilarity = vec![vec![0.0; n]; n];

        // Compute pairwise dissimilarities
        if n >= 10 {
            // Parallel computation for larger matrices
            let results: Vec<_> = (0..n)
                .into_par_iter()
                .flat_map(|i| {
                    ((i + 1)..n)
                        .map(|j| {
                            let d = self.mountford_distance(i, j);
                            (i, j, d)
                        })
                        .collect::<Vec<_>>()
                })
                .collect();

            for (i, j, d) in results {
                self.dissimilarity[i][j] = d;
                self.dissimilarity[j][i] = d;
            }
        } else {
            // Sequential for small matrices
            for i in 0..n {
                for j in (i + 1)..n {
                    let d = self.mountford_distance(i, j);
                    self.dissimilarity[i][j] = d;
                    self.dissimilarity[j][i] = d;
                }
            }
        }
    }

    /// Compute Mountford dissimilarity between two samples
    #[inline]
    fn mountford_distance(&self, i: usize, j: usize) -> f64 {
        let sample_a = &self.presence[i];
        let sample_b = &self.presence[j];

        // Count species in each sample and common species
        let mut a = 0usize; // Species in sample A
        let mut b = 0usize; // Species in sample B
        let mut j_common = 0usize; // Common species

        for k in 0..sample_a.len() {
            let in_a = sample_a[k];
            let in_b = sample_b[k];

            if in_a {
                a += 1;
            }
            if in_b {
                b += 1;
            }
            if in_a && in_b {
                j_common += 1;
            }
        }

        // Mountford index: M = 2J / (2ab - (a+b)J)
        // We return dissimilarity, so we compute 1 - similarity
        // But Mountford is already a similarity index, so we compute 1 - M

        if j_common == 0 {
            return 1.0; // Maximum dissimilarity if no common species
        }

        let a = a as f64;
        let b = b as f64;
        let j = j_common as f64;

        let denominator = 2.0 * a * b - (a + b) * j;

        if denominator <= 0.0 {
            return 0.0; // Perfect similarity
        }

        let similarity = (2.0 * j) / denominator;

        // Clamp to [0, 1] and return dissimilarity
        1.0 - similarity.clamp(0.0, 1.0)
    }

    /// Write dissimilarity matrix to output file
    pub fn output<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = csv::Writer::from_path(path)?;

        // Write header
        let mut header = vec!["".to_string()];
        header.extend(self.samples.clone());
        writer.write_record(&header)?;

        // Write rows
        for (i, sample) in self.samples.iter().enumerate() {
            let mut row = vec![sample.clone()];
            for j in 0..self.samples.len() {
                row.push(format!("{:.6}", self.dissimilarity[i][j]));
            }
            writer.write_record(&row)?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Get the dissimilarity matrix
    pub fn dissimilarity(&self) -> &Vec<Vec<f64>> {
        &self.dissimilarity
    }

    /// Get sample names
    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    /// Get species names
    pub fn species(&self) -> &[String] {
        &self.species
    }

    /// Get number of samples
    pub fn num_samples(&self) -> usize {
        self.samples.len()
    }

    /// Get number of species
    pub fn num_species(&self) -> usize {
        self.species.len()
    }

    /// Create plugin from abundance matrix directly (for testing)
    pub fn from_matrix(samples: Vec<String>, species: Vec<String>, abundance: Vec<Vec<f64>>) -> Self {
        let presence = abundance
            .iter()
            .map(|sample| sample.iter().map(|&v| v > 0.0).collect())
            .collect();

        Self {
            samples,
            species,
            abundance,
            presence,
            dissimilarity: Vec::new(),
        }
    }
}

impl Default for MountfordPlugin {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_csv(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        write!(file, "{}", content).unwrap();
        file
    }

    #[test]
    fn test_identical_samples() {
        // Two identical samples should have dissimilarity 0
        let samples = vec!["A".to_string(), "B".to_string()];
        let species = vec!["sp1".to_string(), "sp2".to_string(), "sp3".to_string()];
        let abundance = vec![
            vec![1.0, 2.0, 3.0], // Sample A
            vec![1.0, 2.0, 3.0], // Sample B (identical)
        ];

        let mut plugin = MountfordPlugin::from_matrix(samples, species, abundance);
        plugin.run();

        assert!((plugin.dissimilarity[0][1] - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_completely_different_samples() {
        // Two samples with no common species should have dissimilarity 1
        let samples = vec!["A".to_string(), "B".to_string()];
        let species = vec![
            "sp1".to_string(),
            "sp2".to_string(),
            "sp3".to_string(),
            "sp4".to_string(),
        ];
        let abundance = vec![
            vec![1.0, 1.0, 0.0, 0.0], // Sample A has sp1, sp2
            vec![0.0, 0.0, 1.0, 1.0], // Sample B has sp3, sp4
        ];

        let mut plugin = MountfordPlugin::from_matrix(samples, species, abundance);
        plugin.run();

        assert!((plugin.dissimilarity[0][1] - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_partial_overlap() {
        // Samples with some common species
        let samples = vec!["A".to_string(), "B".to_string()];
        let species = vec!["sp1".to_string(), "sp2".to_string(), "sp3".to_string()];
        let abundance = vec![
            vec![1.0, 1.0, 0.0], // Sample A has sp1, sp2
            vec![0.0, 1.0, 1.0], // Sample B has sp2, sp3
        ];

        let mut plugin = MountfordPlugin::from_matrix(samples, species, abundance);
        plugin.run();

        // Should have intermediate dissimilarity
        let d = plugin.dissimilarity[0][1];
        assert!(d > 0.0 && d < 1.0);
    }

    #[test]
    fn test_csv_parsing() {
        let csv = "species,Sample1,Sample2,Sample3\nsp1,10,0,5\nsp2,0,20,15\nsp3,5,5,0";
        let file = create_test_csv(csv);

        let mut plugin = MountfordPlugin::new();
        plugin.input(file.path()).unwrap();

        assert_eq!(plugin.num_samples(), 3);
        assert_eq!(plugin.num_species(), 3);
        assert_eq!(plugin.samples(), &["Sample1", "Sample2", "Sample3"]);
    }

    #[test]
    fn test_symmetry() {
        let samples = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let species = vec!["sp1".to_string(), "sp2".to_string(), "sp3".to_string()];
        let abundance = vec![
            vec![1.0, 0.0, 1.0],
            vec![0.0, 1.0, 1.0],
            vec![1.0, 1.0, 0.0],
        ];

        let mut plugin = MountfordPlugin::from_matrix(samples, species, abundance);
        plugin.run();

        // Distance matrix should be symmetric
        for i in 0..3 {
            for j in 0..3 {
                assert!((plugin.dissimilarity[i][j] - plugin.dissimilarity[j][i]).abs() < 0.0001);
            }
        }

        // Diagonal should be 0
        for i in 0..3 {
            assert!((plugin.dissimilarity[i][i] - 0.0).abs() < 0.0001);
        }
    }

    #[test]
    fn test_output_format() {
        let samples = vec!["A".to_string(), "B".to_string()];
        let species = vec!["sp1".to_string()];
        let abundance = vec![vec![1.0], vec![1.0]];

        let mut plugin = MountfordPlugin::from_matrix(samples, species, abundance);
        plugin.run();

        let output_file = NamedTempFile::new().unwrap();
        plugin.output(output_file.path()).unwrap();

        let content = std::fs::read_to_string(output_file.path()).unwrap();
        assert!(content.contains(",A,B"));
        assert!(content.contains("A,"));
        assert!(content.contains("B,"));
    }
}
