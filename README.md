# landmark-data-processing

R scripts for processing geometric morphometric landmark datasets and evaluating density-based missing-landmark estimation.

This repository is associated with the manuscript:

> **Traditional characters and Procrustes-aligned landmark data: a sensitivity analysis in morphological data type weighting for phylogenetic analyses**

A formal citation will be added after publication.

## Citation

If you use this repository, please cite the associated manuscript:

[will be added after publication]

## Overview

This repository contains two R scripts used to process landmark data and evaluate a density-based missing-data estimation procedure:

1. `gm_tps_processing_pipeline.R`
   General TPS landmark-processing pipeline. This script reads and merges TPS files, optionally estimates one-sided missing bilateral landmarks, estimates remaining missing landmarks for GPA alignment, restores the missing-data mask, optionally creates taxon/group mean shapes, and writes a final ready-for-analysis TPS file.

2. `density-based-estimation.R`
   Density-based missing-landmark estimation script. This script calculates pairwise Procrustes distances using shared landmarks, smooths the resulting distance matrix with SMACOF, estimates local morphospace density, applies density-based weights during PCA-based imputation, and compares weighted and unweighted outputs.

## Repository structure

A typical working directory should look like this:

```text
landmark-data-processing/
├── README.md
├── gm_tps_processing_pipeline.R
├── density-based-estimation.R
├── data/
│   └── raw_tps/
│       ├── specimen_set_1.tps
│       └── specimen_set_2.tps
├── metadata/
│   ├── bilateral_pairs.csv
│   └── group_map.csv
└── results/
```

## Requirements

The scripts require R and the following R packages:

```r
library(geomorph)
library(abind)
library(missMDA)
library(readxl)
library(smacof)
library(vegan)
library(FNN)
library(MASS)
```


## Script 1: TPS processing pipeline

### Purpose

`gm_tps_processing_pipeline.R` prepares specimen-level or taxon-level TPS landmark datasets for downstream analysis. It follows the processing order described in the manuscript:

1. read and merge one or more TPS files;
2. optionally estimate one-sided missing bilateral landmarks using `geomorph::bilat.symmetry()`;
3. estimate remaining missing landmarks temporarily to permit global GPA alignment;
4. restore the missing-data mask after alignment;
5. optionally create taxon/group mean shapes using a group map;
6. optionally realign the group-level dataset and restore group-level missing data;
7. write a final ready-for-analysis TPS file.

### Input files

#### TPS files

Place TPS files in:

```text
data/raw_tps/
```

Each TPS file should contain `ID=` specimen identifiers. All input TPS files must have the same number of landmarks and coordinate dimensions.

#### Optional bilateral-pair file

To estimate one-sided missing bilateral landmarks, provide a CSV or Excel file such as:

```text
metadata/bilateral_pairs.csv
```

The default expected columns are:

```csv
landmark,pair
1,2
3,4
5,6
```

Each row identifies a left/right landmark pair. The script removes duplicate reciprocal pairs internally.

#### Optional group map

To collapse multiple specimens into taxon/group means, provide a CSV file such as:

```text
metadata/group_map.csv
```

The file must contain:

```csv
group,specimen_id
Group_A,A1
Group_A,A2
Group_B,B1
Group_B,B2
```

Specimen IDs must match the `ID=` names in the TPS files.

### Configuration

Edit the configuration section at the top of `gm_tps_processing_pipeline.R`.

For example:

```r
input_tps_dir <- file.path("data", "raw_tps")
input_tps_files <- NULL
output_dir <- file.path("results", "gm_tps_processing_pipeline")

bilateral_pairs_file <- file.path("metadata", "bilateral_pairs.csv")
group_map_file <- file.path("metadata", "group_map.csv")
```

Use `input_tps_files <- NULL` to read all `.tps` files in `input_tps_dir`. Alternatively, provide an explicit file list:

```r
input_tps_files <- c(
  file.path("data", "raw_tps", "specimen_set_1.tps"),
  file.path("data", "raw_tps", "specimen_set_2.tps")
)
```

### Outputs

The script creates the following output structure:

```text
results/gm_tps_processing_pipeline/
├── intermediate/
├── group_means/
├── ready_for_analysis/
└── tables/
```

`intermediate/` contains stepwise TPS outputs, including the merged raw dataset, the post-bilateral-symmetry dataset, the globally aligned dataset with missing values restored, and the optional group/taxon mean dataset.

`group_means/` contains individual TPS files for each group/taxon mean when a group map is supplied.

`ready_for_analysis/` contains the final TPS files for downstream analyses. The script writes both a timestamped file and a stable `ready_for_analysis_latest.tps` file.

`tables/` contains CSV summaries documenting missingness, bilateral-symmetry filling, imputation settings, group composition, output paths, and session information.

## Script 2: Density-based missing-landmark estimation

### Purpose

`density-based-estimation.R` evaluates whether uneven sampling density in morphospace affects missing landmark estimation. Be careful in how this is applied as it may obscure real biological signal.

The script:

1. reads a TPS landmark dataset;
2. calculates pairwise Procrustes distances using only mutually observed landmarks;
3. corrects distances by the proportion of shared landmarks;
4. smooths the corrected distance matrix using SMACOF;
5. quantifies clustering using Hopkins’ statistic across repeated random subsamples;
6. estimates local morphospace radius using nearest neighbors;
7. defines density weights proportional to local radius and normalized to mean one;
8. performs weighted and unweighted PCA-based missing-landmark estimation;
9. compares weighted and unweighted outputs after GPA and NMDS.

### Configuration

Edit the input and output paths near the top of `density-based-estimation.R`, or pass them from the command line depending on the version of the script.

Example:

```r
input_tps <- file.path("results", "gm_tps_processing_pipeline", "ready_for_analysis", "ready_for_analysis_latest.tps")
output_dir <- file.path("results", "density_based_estimation")
```

### Outputs

The density-based estimation script writes diagnostic tables and figures, including:

```text
input_missingness_summary.csv
pairwise_procrustes_distances_raw.csv
pairwise_procrustes_distances_corrected.csv
pairwise_shared_landmarks.csv
smacof_coordinates.csv
hopkins_statistic.csv
hopkins_statistic_replicates.csv
density_weights.csv
density_weight_summary.csv
pca_imputation_summary.csv
nmds_weighted_unweighted_coordinates.csv
density_weighting_diagnostics_summary.csv
```

It also writes complete weighted and unweighted TPS outputs for comparison.

## Notes on missing data

Both scripts temporarily estimate missing landmarks for specific computational steps. The TPS processing script restores the appropriate missing-data mask after GPA alignment and after group-level processing. In the density-based estimation script, weighted and unweighted complete configurations are intentionally written so the effect of weighting can be compared directly.

## Notes on bilateral symmetry

The TPS processing script uses `geomorph::bilat.symmetry()` when a bilateral-pair file is supplied. The script fills only landmarks where one member of a bilateral pair is missing and the other is observed. Observed landmarks are not replaced, and pairs where both landmarks are missing remain missing.
