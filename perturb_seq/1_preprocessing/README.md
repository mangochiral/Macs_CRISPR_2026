# 1. Preprocessing

This directory contains scripts and instructions for preprocessing the perturb-seq data prior to downstream analysis. Proper preprocessing is essential to ensure high-quality, reliable results in subsequent steps.

## Overview

The preprocessing workflow includes:
- Ingestion of raw data and spliting it into CRISPR and Gene Expression anndata objects
- Quality control checks on raw data
- Filtering of low-quality cells and genes
- Normalization of gene expression counts
- Computation total guide UMI counts for each cell.
- Conversion of data into analysis-ready formats

## Folder Structure

```
.
├── 1_preprocessing/
│   ├── preprocess_macs.py
│   ├── preprocess_adata.py                # importing modules
│   ├── basic_preprocessing.py       # for parallel run
```

> [`preprocess_adata.py`](https://github.com/mangochiral/PerturbSeq_Analysis_pipeline](https://github.com/mangochiral/PerturbSeq_Analysis_pipeline/blob/main/src/1_preprocessing/preprocess_adata.py) is sourced from the PerturbSeq Analysis Pipeline repository.

## Input Data
The experiment_info.csv has sample info about condition, lane, donor, treatment, etc
sample_filter_feature_matrix.h5d is cellranger output
- file1: sample_filter_feature_matrix.h5d
- file2: experiment_info.csv
  
## Output files:
- file1: `<sample_name>_gex_preprocessed.h5ad`
- file2: `<sample_name>_crispr_preprocessed.h5ad`

Initial QC was done to remove low quality cells i.e, high mitochondrial percent (> 20%, decent cut for macrophage cells), and cells have less than 200 genes.

## Running the Pipeline

1. Adjust paths and parameters as needed in the scripts.
2. From this directory, run the main preprocessing script:
   scripts 
   bash script for slurm job
   ```
    #!/bin/bash
    #SBATCH --job-name=preprocessing
	#SBATCH --time=1:00:00
	#SBATCH --mem=100G
	#SBATCH --cpus-per-task=16         # <-- 16 CPUs per job
	#SBATCH --output=logs/preprocess_%A_%a.out
	#SBATCH --error=logs/preprocess_%A_%a.err


	# Use single-threaded BLAS per process (safer with multiprocessing)
	export OMP_NUM_THREADS=1

	# Run the guide assignment
	python3 basic_preprocessing_macs.py  --datadir <PATH TO CELLRANGER OUTPUT> \
  	--nprocs "${SLURM_CPUS_PER_TASK}"   --output_dir <PATH TO PROCESSED DATA OUTPUT> --exp_info <EXPERIMENT METADATA> --exp 'crispr'
	echo "Completed!"

   ````