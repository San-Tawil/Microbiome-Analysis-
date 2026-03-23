# Microbiome Analysis Pipeline

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A complete, reproducible microbiome analysis pipeline for MetaPhlAn-derived relative abundance tables. Performs taxonomic profiling, alpha & beta diversity analysis, ordination, and permutational statistical testing.

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation & Environment Setup](#installation--environment-setup)
- [Input Data](#input-data)
- [Step-by-Step Usage Guide](#step-by-step-usage-guide)
- [Output Structure](#output-structure)
- [Methods](#methods)
- [References](#references)
- [Citation](#citation)

---

## Overview

This pipeline takes pre-computed relative abundance tables (e.g., from MetaPhlAn [[1]](#ref-metaphlan) or QIIME 2 [[2]](#ref-qiime2)) at multiple taxonomic ranks and a sample metadata file, then performs:

1. **Taxonomic profiling** — horizontal stacked barplots grouped by condition and per sample.
2. **Alpha diversity** — Shannon, Simpson, and Inverse Simpson indices with statistical testing.
3. **Beta diversity** — Bray–Curtis, Jaccard, and Aitchison distance matrices with PCoA and NMDS ordination.
4. **Statistical testing** — PERMANOVA and PERMDISP with 999 permutations.

All outputs (high-resolution plots, CSV tables, distance matrices, ordination coordinates, and test results) are saved to a structured `results/` directory. The script uses the `Agg` backend for Matplotlib, enabling execution on headless servers and HPC clusters.

---

## Requirements

### Python Version

- Python ≥ 3.8

### Dependencies

| Package | Min. Version | Purpose |
|---------|-------------|---------|
| `numpy` [[3]](#ref-numpy) | 1.21 | Numerical computation |
| `pandas` [[4]](#ref-pandas) | 1.3 | Data manipulation |
| `matplotlib` [[5]](#ref-matplotlib) | 3.4 | Plotting (Agg backend for headless) |
| `seaborn` [[6]](#ref-seaborn) | 0.11 | Statistical visualization |
| `scipy` [[7]](#ref-scipy) | 1.7 | Distance metrics, statistical tests |
| `scikit-learn` [[8]](#ref-sklearn) | 1.0 | NMDS via multidimensional scaling |

---

## Installation & Environment Setup

### Option A: pip (virtualenv)

```bash
# 1. Create and activate a virtual environment
python -m venv microbiome_env

# Linux / macOS:
source microbiome_env/bin/activate
# Windows:
microbiome_env\Scripts\activate

# 2. Install dependencies
pip install numpy pandas matplotlib seaborn scipy scikit-learn
```

### Option B: Conda

```bash
# 1. Create the environment
conda create -n microbiome python=3.10 -y

# 2. Activate
conda activate microbiome

# 3. Install dependencies
conda install numpy pandas matplotlib seaborn scipy scikit-learn -y
```

---

## Input Data

Place all input files in the **same directory** as `microbiome_analysis.py`.

### 1. Relative Abundance Tables

Four tab-separated files, one per taxonomic level:

| Filename | Taxonomic Level |
|----------|----------------|
| `8_rel-phyla-table.tsv` | Phylum |
| `8_rel-family-table.tsv` | Family |
| `8_rel-genus-table.tsv` | Genus |
| `8_rel-species-table.tsv` | Species |

**File format:**

```
# Constructed from biom file
#OTU ID        Sample1  Sample2  Sample3  ...
d__Bacteria;k__Bacteria;p__Firmicutes     0.45   0.32   0.51   ...
d__Bacteria;k__Bacteria;p__Bacteroidota   0.30   0.41   0.28   ...
```

- **Row 1**: Comment line (skipped automatically).
- **Row 2**: Header — `#OTU ID` followed by sample identifiers.
- **Rows 3+**: Semicolon-delimited taxonomy string, then relative abundance values (proportions ∈ [0, 1]).

### 2. Metadata Table

A tab-separated file named `metadata.tsv`:

```
#SampleID   State
FMF-C1      Control
FMF-P1      FMF
...
```

- **First column**: Sample IDs (must *exactly* match column headers in abundance tables).
- **`State` column**: Grouping variable for all comparisons (e.g., `FMF` vs `Control`).

> **Note:** The pipeline validates sample ID overlap between tables and metadata and warns about any mismatches.

---

## Step-by-Step Usage Guide

### Step 1: Prepare Your Data

1. Place the four relative abundance `.tsv` files and `metadata.tsv` in a single directory.
2. Verify that sample IDs in the abundance table headers match the `#SampleID` column in `metadata.tsv`.
3. Ensure the `State` column in your metadata contains the grouping labels you want to compare (e.g., `FMF` and `Control`).

### Step 2: Set Up the Environment

```bash
# Create and activate a virtual environment (if not done already)
python -m venv microbiome_env

# Windows:
microbiome_env\Scripts\activate
# Linux/macOS:
source microbiome_env/bin/activate

# Install all required packages
pip install numpy pandas matplotlib seaborn scipy scikit-learn
```

### Step 3: Configure the Script (Optional)

Open `microbiome_analysis.py` and edit the configuration section near the top of the file:

```python
# --- Configurable parameters ---
GROUP_COL    = 'State'      # Metadata column for grouping
TOP_N_VALUES = [10, 15]     # Top N taxa to display in plots
RANDOM_SEED  = 42           # Seed for NMDS reproducibility
```

If your abundance files have different names, update the `TABLE_FILES` dictionary:

```python
TABLE_FILES = {
    'phylum':  os.path.join(SCRIPT_DIR, 'your_phylum_file.tsv'),
    'family':  os.path.join(SCRIPT_DIR, 'your_family_file.tsv'),
    'genus':   os.path.join(SCRIPT_DIR, 'your_genus_file.tsv'),
    'species': os.path.join(SCRIPT_DIR, 'your_species_file.tsv'),
}
```

### Step 4: Run the Pipeline

```bash
cd /path/to/data/directory
python microbiome_analysis.py
```

The script will log each step to the console:

```
01:34:10 | INFO | ============================================
01:34:10 | INFO | MICROBIOME ANALYSIS PIPELINE
01:34:10 | INFO | ============================================
01:34:10 | INFO | Metadata: 29 samples, columns=['State']
01:34:10 | INFO | [phylum] 29 samples × 27 taxa
01:34:11 | INFO | STEP 1: Taxonomic profiling
01:34:15 | INFO | STEP 2: Alpha diversity
01:34:17 | INFO | STEP 3: Beta diversity & statistics
01:34:21 | INFO | ============================================
01:34:21 | INFO | COMPLETE — results in: .../results
01:34:21 | INFO | ============================================
```

### Step 5: Examine the Results

All outputs are saved in the `results/` subdirectory. Open PNG files with any image viewer and CSV files with Excel, R, or pandas:

```bash
# List all generated files
# Linux/macOS:
find results/ -type f | sort

# Windows:
dir /s /b results\
```

### Step 6: Run on a Remote Server / HPC (Optional)

The script uses `matplotlib.use('Agg')` — no display server or X11 forwarding is needed:

```bash
ssh user@server
cd /path/to/data
python microbiome_analysis.py   # No X11 needed
```

---

## Output Structure

```
results/
├── taxa_barplots/
│   ├── {level}_top{N}_barplot.png              # Grouped barplot (Patient vs Control)
│   ├── {level}_top{N}_per_sample_barplot.png   # Per-sample barplot
│   ├── {level}_top{N}_group_means.csv          # Group-level mean abundances (%)
│   └── {level}_top{N}_per_sample_data.csv      # Per-sample abundances (%)
│
├── alpha_diversity/
│   ├── plots/
│   │   └── {level}_alpha_boxplots.png          # Boxplots per diversity metric
│   ├── values/
│   │   └── {level}_alpha_diversity.csv         # Per-sample diversity values
│   └── stats/
│       └── alpha_diversity_stats.csv           # Kruskal–Wallis & Mann–Whitney U
│
└── beta_diversity/
    ├── distance_matrices/
    │   └── {metric}_distance_matrix.csv        # Full pairwise distance matrix
    ├── pcoa/
    │   ├── PCoA_{metric}_top{N}.png            # PCoA scatter with 95% ellipses
    │   └── PCoA_{metric}_coords.csv            # PCoA coordinates
    ├── nmds/
    │   ├── NMDS_{metric}_top{N}.png            # NMDS scatter with 95% ellipses
    │   └── NMDS_{metric}_coords.csv            # NMDS coordinates
    ├── boxplots/
    │   ├── {metric}_boxplot.png                # Within vs between group distances
    │   └── {metric}_boxplot_data.csv
    ├── permanova/
    │   └── PERMANOVA_{metric}.csv              # Per-metric PERMANOVA result
    ├── permdisp/
    │   └── PERMDISP_{metric}.csv               # Per-metric PERMDISP result
    └── beta_stats_summary.csv                  # Consolidated statistics
```

Where:
- `{level}` ∈ {`phylum`, `family`, `genus`, `species`}
- `{metric}` ∈ {`Bray_Curtis`, `Jaccard`, `Aitchison`}
- `{N}` ∈ {`10`, `15`}

---

## Methods

### Taxonomic Profiling

Relative abundances are visualized as horizontal stacked barplots. "Unassigned" and "Other" categories are excluded from visualization; the remaining taxa are renormalized so that each sample or group sums to 100%. Taxa are ranked by global mean relative abundance and the top *N* are displayed with a consistent colour mapping across all plots.

### Alpha Diversity

Three diversity indices are computed per sample from relative abundance vectors:

- **Shannon index** [[9]](#ref-shannon): $H' = -\sum_{i=1}^{S} p_i \ln(p_i)$  
  where *S* is the number of taxa and *p*ᵢ is the relative abundance of taxon *i*.
- **Simpson index** [[10]](#ref-simpson): $D = \sum_{i=1}^{S} p_i^2$
- **Inverse Simpson**: $D^{-1} = 1 / D$

Between-group comparisons use the Kruskal–Wallis *H*-test [[11]](#ref-kruskal) and pairwise Mann–Whitney *U*-test [[12]](#ref-mann).

### Beta Diversity

Three distance metrics are computed from the genus-level abundance table:

| Metric | Input | Formula |
|--------|-------|---------|
| **Bray–Curtis** [[13]](#ref-bray) | Relative abundance | $BC_{jk} = 1 - \frac{2 \sum \min(x_{ij}, x_{ik})}{\sum x_{ij} + \sum x_{ik}}$ |
| **Jaccard** [[14]](#ref-jaccard) | Presence/absence | $J = 1 - \frac{|A \cap B|}{|A \cup B|}$ |
| **Aitchison** [[15]](#ref-aitchison) | CLR-transformed | Euclidean distance on CLR values (pseudocount 10⁻⁶) |

The centred log-ratio (CLR) transformation:

$$\text{clr}(\mathbf{x})_i = \ln(x_i) - \frac{1}{S}\sum_{j=1}^{S}\ln(x_j)$$

A pseudocount of 10⁻⁶ is added before log-transformation to handle zero values.

### Ordination

- **PCoA** (Principal Coordinates Analysis) [[16]](#ref-gower): classical multidimensional scaling via eigendecomposition of the doubly-centred distance matrix.
- **NMDS** (Non-metric Multidimensional Scaling) [[17]](#ref-kruskal-mds): performed via `sklearn.manifold.MDS` with `random_state = 42` for reproducibility.

95% confidence ellipses are drawn per group on all ordination plots.

### Statistical Tests

- **PERMANOVA** [[18]](#ref-anderson2001): Permutational Multivariate Analysis of Variance (999 permutations). Tests whether group centroids differ significantly in multivariate space:

$$F = \frac{SS_A / (k-1)}{SS_W / (n-k)}$$

  where *SS*_A and *SS*_W are the among- and within-group sums of squared distances, *k* is the number of groups, and *n* the total number of samples.

- **PERMDISP** [[19]](#ref-anderson2006): tests homogeneity of multivariate dispersions (i.e., whether the spread of samples around each group centroid is equal) using 999 permutations.

---

## References

### Methods

<a id="ref-shannon"></a>[9] Shannon, C. E. (1948). A mathematical theory of communication. *The Bell System Technical Journal*, 27(3), 379–423. [doi:10.1002/j.1538-7305.1948.tb01338.x](https://doi.org/10.1002/j.1538-7305.1948.tb01338.x)

<a id="ref-simpson"></a>[10] Simpson, E. H. (1949). Measurement of diversity. *Nature*, 163(4148), 688. [doi:10.1038/163688a0](https://doi.org/10.1038/163688a0)

<a id="ref-kruskal"></a>[11] Kruskal, W. H. & Wallis, W. A. (1952). Use of ranks in one-criterion variance analysis. *J. Am. Stat. Assoc.*, 47(260), 583–621. [doi:10.1080/01621459.1952.10483441](https://doi.org/10.1080/01621459.1952.10483441)

<a id="ref-mann"></a>[12] Mann, H. B. & Whitney, D. R. (1947). On a test of whether one of two random variables is stochastically larger than the other. *Ann. Math. Stat.*, 18(1), 50–60. [doi:10.1214/aoms/1177730491](https://doi.org/10.1214/aoms/1177730491)

<a id="ref-bray"></a>[13] Bray, J. R. & Curtis, J. T. (1957). An ordination of the upland forest communities of southern Wisconsin. *Ecological Monographs*, 27(4), 325–349. [doi:10.2307/1942268](https://doi.org/10.2307/1942268)

<a id="ref-jaccard"></a>[14] Jaccard, P. (1912). The distribution of the flora in the alpine zone. *New Phytologist*, 11(2), 37–50. [doi:10.1111/j.1469-8137.1912.tb05611.x](https://doi.org/10.1111/j.1469-8137.1912.tb05611.x)

<a id="ref-aitchison"></a>[15] Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*. Chapman and Hall. [doi:10.1007/978-94-009-4109-0](https://doi.org/10.1007/978-94-009-4109-0)

<a id="ref-gower"></a>[16] Gower, J. C. (1966). Some distance properties of latent root and vector methods used in multivariate analysis. *Biometrika*, 53(3–4), 325–338. [doi:10.1093/biomet/53.3-4.325](https://doi.org/10.1093/biomet/53.3-4.325)

<a id="ref-kruskal-mds"></a>[17] Kruskal, J. B. (1964). Multidimensional scaling by optimizing goodness of fit to a nonmetric hypothesis. *Psychometrika*, 29(1), 1–27. [doi:10.1007/BF02289565](https://doi.org/10.1007/BF02289565)

<a id="ref-anderson2001"></a>[18] Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology*, 26(1), 32–46. [doi:10.1111/j.1442-9993.2001.01070.pp.x](https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x)

<a id="ref-anderson2006"></a>[19] Anderson, M. J. (2006). Distance-based tests for homogeneity of multivariate dispersions. *Biometrics*, 62(1), 245–253. [doi:10.1111/j.1541-0420.2005.00440.x](https://doi.org/10.1111/j.1541-0420.2005.00440.x)

### Bioinformatics Tools

<a id="ref-metaphlan"></a>[1] Blanco-Míguez, A., et al. (2023). Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nature Biotechnology*, 41, 1633–1644. [doi:10.1038/s41587-023-01688-w](https://doi.org/10.1038/s41587-023-01688-w)

<a id="ref-qiime2"></a>[2] Bolyen, E., et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology*, 37, 852–857. [doi:10.1038/s41587-019-0209-9](https://doi.org/10.1038/s41587-019-0209-9)

### Software

<a id="ref-numpy"></a>[3] Harris, C. R., et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362. [doi:10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)

<a id="ref-pandas"></a>[4] McKinney, W. (2010). Data structures for statistical computing in Python. *Proc. 9th Python in Science Conf.*, 56–61. [doi:10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)

<a id="ref-matplotlib"></a>[5] Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Comput. Sci. Eng.*, 9(3), 90–95. [doi:10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55)

<a id="ref-seaborn"></a>[6] Waskom, M. L. (2021). seaborn: statistical data visualization. *JOSS*, 6(60), 3021. [doi:10.21105/joss.03021](https://doi.org/10.21105/joss.03021)

<a id="ref-scipy"></a>[7] Virtanen, P., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272. [doi:10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)

<a id="ref-sklearn"></a>[8] Pedregosa, F., et al. (2011). Scikit-learn: Machine Learning in Python. *JMLR*, 12, 2825–2830. [jmlr.org/papers/v12/pedregosa11a](https://jmlr.org/papers/v12/pedregosa11a.html)

---

## Citation

If you use this pipeline in your research, please cite the relevant methods and software listed in the [References](#references) section above.

---
