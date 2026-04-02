# MACS KO Pseudobulk DE Pipeline

Pseudobulking was done by separating the conditions. This process helps in simplifying the DEseq2 design for model fitting.
---

## Pipeline Overview
Source code: [src](https://github.com/mangochiral/PerturbSeq_Analysis_pipeline/tree/main/src/3_pseudobulk)

```
(Optional) Thresholding_doublet_in_multiguide.ipynb

Step 1)  Keep_singlets_prep_adata.ipynb
Step 2)  pseudobulk_by_lane.py
Step 3)  prep_DE_MACS.ipynb
Step 4)  Deseq2_pseudobulk.py

Downstream visualization
  ├── DE_results_plot.ipynb
  └── Regulator_heatmap.ipynb
```

---

## (Optional) Doublet Removal in Multi-guide Cells (`Thresholding_doublet_in_multiguide.ipynb`)

If multi-guide–assigned cells are used, this notebook should be run **prior to Step 1** to identify and remove likely doublets from the multi-guide pool.

### What it does
**Step 1 - `Keep_singlets_prep_adata.ipynb`: Filter to singlets**: Remove cells with no guide, multiple guides, or ambiguous assignments. Keep only single-sgRNA and single-NTC cells.

### Inputs / Outputs

| | File | Description |
|---|---|---|
| **In** | `<sample_name>_gex_guide.h5ad` |
| **Out** | `<sample_name>_singlets.h5ad` |


**Step 2 — Per-lane pseudobulk**: Aggregate filtered cells into pseudobulk profiles (sum of raw counts) grouped by condition-level covariates.

```bash
#!/bin/bash
#SBATCH --job-name=pseudo_bulk
#SBATCH --array=0-7
#SBATCH --time=2:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/pseudo_%A_%a.out
#SBATCH --error=logs/pseudo_%A_%a.err


# Safer with numpy/BLAS + multiprocessing
export OMP_NUM_THREADS=1

PROCESSED_DIR=<PATH TO PROCESSED DATA DIRECTORY>
EXPMETA=<PATH TO EXPERIMENT METADATA CSV>


python3 pseudobulk_by_lane.py \
  --processed_dir "${PROCESSED_DIR}" \
  --expmeta "${EXPMETA}" \
  --nprocs "${SLURM_CPUS_PER_TASK}" \
  --task_id "${SLURM_ARRAY_TASK_ID}"\
  --sgrna_col "assigned_guide_id"

```

### Inputs / Outputs

| | File | Description |
|---|---|---|
| **In** | `<sample_name>_singlets.h5ad`  
| **Out** | `<sample_name>_DE_pseudobulk.h5ad`


**Step 3 — Merge, QC & prep for DE**: Concatenate per-lane pseudobulk files, apply quality filters (minimum cells, total-count outliers, replicate thresholds), select features, and write the final DE-ready h5ad.

### Inputs / Outputs

| | File | Description |
|---|---|---|
| **In** | `<sample_name>_DE_pseudobulk.h5ad` |
| **Out** | `merged_DE_pseudobulk.h5ad` |
| **Out** | `<condition>_DE_pseudobulk_for_test.h5ad` | 
| **Out** | `DE_feature_selection_vars.csv` | 
---

## Step 4 — DESeq2 on Pseudobulk ([Deseq2_pseudobulk.py](https://github.com/mangochiral/PerturbSeq_Analysis_pipeline/tree/main/src/4_DGE_analysis))

### What it does

1. Loads the DE-ready pseudobulk h5ad (`*_DE_pseudobulk_for_test.h5ad`) produced by Step 3.
2. Subsets to highly variable genes identified during feature selection (`DE_feature_selection_vars.csv`).
3. For each condition, fits a negative-binomial GLM via `pertpy.tl.PyDESeq2` with the design `~ log10_n_cells + target_gene`, regressing out cell-count as a confounder.
5. For each target gene in a chunk, performs a Wald contrast (target vs NTC) followed by FDR correction.
6. Chunk-level CSVs are saved to `DE_outputs/` and concatenated into a single results table.
7. Results are pivoted into an AnnData where observations are contrasts, variables are genes, and layers store `baseMean`, `log_fc`, `lfcSE`, `p_value`, and `adj_p_value`. Written as `DE_<condition>_anndata.h5ad`.

### Usage

```bash
#!/bin/bash
#SBATCH --job-name=guide_assign
#SBATCH --array=0
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/deseq_macs_%A_%a.out
#SBATCH --error=logs/deseq_macs_%A_%a.err

export OMP_NUM_THREADS=1

PROCESSED_DIR=<PATH TO PROCESSED DATA DIRECTORY>


# Derive sample and condition from array index

python3 Deseq2_pseudobulk.py \
  --processed_dir "${PROCESSED_DIR}"  \
  --nprocs "${SLURM_CPUS_PER_TASK}" \
 --task_id "${SLURM_ARRAY_TASK_ID}"
echo "Done"


```

### Inputs / Outputs

| | File | Description |
|---|---|---|
| **In** | `<condition>_DE_pseudobulk_for_test.h5ad` | DE-ready pseudobulk (Step 3 output) |
| **Out** | `DE_<condition>_anndata.h5ad` | AnnData with DE statistics in layers |

---

## Downstream Analysis & Visualization

### DE Results Summary (`DE_results_plot.ipynb`)

Explores the DE results produced in Step 4.

1. Loads per-condition DE AnnData objects and computes a z-score layer (`log_fc / lfcSE`).
2. **Per-perturbation DE summary**: Counts the number of up-regulated, down-regulated, and total DE genes per perturbation. Plots ranked scatter plots of total DE genes per perturbation, annotating the top 10 hits.
3. **Heatmap on genes of interest**: Extracts DE results for a curated gene list (cytokines/chemokines, M1/M2 markers, ISGs, phagocytosis markers) and plots z-score cluster-heatmaps with significance dots (FDR < 0.1) and category-colored row annotations, separately for each condition.
4. **Singlet vs multiplet concordance** (if multi-guide cells were retained): Runs Wilcoxon tests on single-cell expression of genes of interest in multiplet-assigned cells (target vs other) across LPS lanes, then correlates the multiplet log-fold-changes against the singlet pseudobulk log-fold-changes gene by gene (Pearson R scatter plots). Assesses whether multiplet-derived effects agree with singlet-based DE.
5. **Regulator identification**: Identifies perturbations that significantly alter any gene of interest (excluding on-target effects), then retrieves their genome-wide DE profiles and plots the number of significant effects per regulator as bar charts.

### Regulator Heatmap (`Regulator_heatmap.ipynb`)

Focused visualization of perturbation effects on the curated gene panel.

1. Loads per-condition DE AnnData and computes z-scores.
2. Defines a gene list organized by functional category (cytokines/chemokines, M1 markers, M2 markers, ISGs, phagocytosis receptors, etc.) and applies manual edits (gene additions/removals).
3. Extracts long-format DE results for the gene panel and plots publication-quality cluster-heatmaps (z-score coloring, `bwr` palette, category row colors) using `seaborn.clustermap`. Perturbations with no significant effects are filtered out. Self-target effects (gene = perturbation) are masked.
4. Generates separate heatmaps for LPS and Ctrl conditions with configurable z-score bounds, dot overlays, and row/column ordering.

