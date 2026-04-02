# Guide Efficiency Analysis

Analysis notebooks for quality control and guide efficiency assessment of a CRISPR Perturb-seq screen in macrophages (Ctrl vs LPS-stimulated conditions).

## Guide Assignment
Gaussian-Poisson mixed model was used to assign guide to cells. The cut off umi counts cut limits for a individual guides are recorded in csv file and for individual cells umi counts for guide assigned are recorded, (for multi guide cells, the guide with highest umi counts is recorded).
> guide_assignment: [guide_assignment.py](https://github.com/mangochiral/PerturbSeq_Analysis_pipeline/blob/main/src/2_guide-assignment/guide_assignment_parallel.py)

Input files : `<sample_name>_crispr_preprocessed.h5ad`
Output files: 	- file1: `{sample_name}_gex_guide.h5ad`
		- file2: `guide_assignment.csv`
		- file1: `guide_threshold.csv`
		- file1: `{sample_name_processed_guide.csv`


## Quality control: To check for guide efficiency, expression levels of perturbed genes were used to test statistics
> [guide qc](https://github.com/mangochiral/PerturbSeq_Analysis_pipeline/blob/main/src/2_guide-assignment/qc_stats_heavy_load.py)

Input files :`{sample_name}_gex_guide.h5ad`
Output files: `{sample_name}_guide_count_info.csv`

## Interactive visualization of Guide efficiency

## Data

All notebooks expect preprocessed `.h5ad` (AnnData) and `.csv` files organized under a shared project directory, with each sequencing lane stored in its own subfolder (`Mac_*`). Key input files per lane include:

- `*_gex_guide.h5ad` — gene expression + guide assignment per cell
- `*_crispr_preprocessed.h5ad` — CRISPR library UMI counts
- `*_guide_count_info.csv` — aggregate guide-level count statistics
- `*_cells_based_guide_efficiency.csv` — cell-level knockdown calls
- `experiment_info.csv` — lane/condition metadata

## Notebooks

### 1. `Guides_per_cell_distribution_plots.ipynb`

Visualizes how many sgRNA guides are detected per cell across conditions.

**What it does:**
- Reads guide assignment tables from each lane and splits by condition (Ctrl / LPS).
- Counts the number of guides assigned to each cell barcode.
- Generates a grouped bar chart (via `plotnine`) comparing the per-cell guide count distribution between Ctrl and LPS.

**Key output:** Bar plot of guide-per-cell distribution, saved to the `plots/` directory.

---

### 2. `Guide_type_distribuition_macs_ko.ipynb`

Characterizes the sgRNA assignment landscape and provides a broad QC overview of the screen.

**What it does:**
- Classifies every cell into an sgRNA category: *single sgRNA*, *single NTC sgRNA*, *multiple sgRNAs*, or *no sgRNA*.
- Plots stacked bar charts of sgRNA type composition per condition.
- Merges all lanes into a single AnnData object and generates dot plots of macrophage marker genes (e.g., CD14, CD68, IL1B, CCL2, CD163) across conditions.
- Computes per-condition summary statistics: mean mRNA UMIs/cell, mean genes/cell, mean % mitochondrial reads, mean guide UMI counts, and mean cells per guide.
- Produces a multi-panel bar chart summarizing these QC metrics.

**Key outputs:** sgRNA-type bar charts, marker gene dot plots, and a QC summary figure.

---

### 3. `QC_stats_MACS_with_cell_level_guide_efficiency.ipynb`

Quantifies guide-level and cell-level knockdown efficiency with statistical testing.

**What it does:**
- Aggregates guide count statistics across lanes (separately for Ctrl and LPS).
- For each targeting guide, computes a Welch's t-test of guide expression vs. NTC expression, derives log₂ fold-change, and applies Benjamini–Hochberg FDR correction.
- Calculates cell-level knockdown rates: the fraction of cells carrying each guide where expression falls below the 5th percentile of the NTC distribution.
- Merges bulk and cell-level efficiency metrics and exports per-condition CSV files.
- Ranks target genes by baseline (NTC) expression and bins them to assess whether knockdown detection depends on expression level.
- Generates scatter plots of guide effect size vs. rank, cell-level knockdown fraction, and binned knockdown efficiency curves.

**Key outputs:** Guide efficiency CSVs (`Macs_Ctrl_cell_level_guide_efficiency.csv`, `Macs_LPS_cell_level_guide_efficiency.csv`), scatter plots, and rank-binned efficiency figures.

## Usage

1. Update the `path_macs` / `lane_path` variables at the top of each notebook to point to your local project directory.
2. Run the notebooks in order (1 → 2 → 3) or independently — they share the same input data but do not depend on each other's outputs.
3. Plots are saved to `plots/` and QC tables to `qc_stats/` within the project directory.
