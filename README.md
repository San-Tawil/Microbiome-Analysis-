# Microbiome Analysis Pipeline

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A complete, reproducible microbiome analysis pipeline for MetaPhlAn-derived relative abundance tables. Performs taxonomic profiling, alpha & beta diversity analysis, ordination, and permutational statistical testing.

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Data](#input-data)
- [Usage](#usage)
- [Output Structure](#output-structure)
- [Methods](#methods)
- [References](#references)
- [Citation](#citation)

---

## Overview

This pipeline takes pre-computed relative abundance tables (e.g., from MetaPhlAn or QIIME 2) at multiple taxonomic ranks and a sample metadata file, then performs:

1. **Taxonomic profiling** — horizontal stacked barplots grouped by condition and per sample
2. **Alpha diversity** — Shannon, Simpson, and Inverse Simpson indices with statistical testing
3. **Beta diversity** — Bray–Curtis, Jaccard, and Aitchison distance matrices with PCoA/NMDS ordination
4. **Statistical testing** — PERMANOVA and PERMDISP with 999 permutations

All outputs (plots, CSV tables, distance matrices, ordination coordinates, and test results) are saved to a structured `results/` directory.

---

## Requirements

### Python Version

- Python ≥ 3.8

### Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `numpy` | ≥ 1.21 | Numerical computation |
| `pandas` | ≥ 1.3 | Data manipulation |
| `matplotlib` | ≥ 3.4 | Plotting (Agg backend for headless) |
| `seaborn` | ≥ 0.11 | Statistical visualization |
| `scipy` | ≥ 1.7 | Distance metrics, statistical tests |
| `scikit-learn` | ≥ 1.0 | NMDS via multidimensional scaling |

### Environment Setup

```bash
# Option 1: pip
pip install numpy pandas matplotlib seaborn scipy scikit-learn

# Option 2: conda
conda create -n microbiome python=3.10
conda activate microbiome
conda install numpy pandas matplotlib seaborn scipy scikit-learn
```

---

## Input Data

The pipeline expects **three types of files** in the same directory as the script:

### 1. Relative Abundance Tables

Tab-separated files with a BIOM-style header. One file per taxonomic level:

| Filename | Taxonomic Level |
|----------|----------------|
| `8_rel-phyla-table.tsv` | Phylum |
| `8_rel-family-table.tsv` | Family |
| `8_rel-genus-table.tsv` | Genus |
| `8_rel-species-table.tsv` | Species |

**Format:**

```
# Constructed from biom file
#OTU ID	Sample1	Sample2	Sample3	...
d__Bacteria;k__Bacteria;p__Firmicutes	0.45	0.32	0.51	...
d__Bacteria;k__Bacteria;p__Bacteroidota	0.30	0.41	0.28	...
```

- **Row 1**: Comment line (skipped automatically)
- **Row 2**: Header with `#OTU ID` followed by sample identifiers
- **Rows 3+**: Semicolon-delimited taxonomy string, then relative abundance values per sample
- Values should be proportions (0–1) summing to ~1 per sample

### 2. Metadata Table

A tab-separated file (`metadata.tsv`) mapping sample IDs to experimental groups:

```
#SampleID	State
FMF-C1	Control
FMF-P1	FMF
...
```

- **First column**: Sample IDs (must match column headers in abundance tables)
- **`State` column**: Grouping variable used for all comparisons (e.g., `FMF` vs `Control`)

> **Note:** Sample IDs in the metadata must exactly match those in the abundance table headers. The pipeline validates this overlap and warns about mismatches.

---

## Usage

### Basic Run

```bash
cd /path/to/data/directory
python microbiome_analysis.py
```

The script automatically discovers input files in its own directory and writes results to `results/`.

### Configuration

Edit the following constants at the top of `microbiome_analysis.py`:

```python
GROUP_COL = 'State'          # Metadata column for group comparison
TOP_N_VALUES = [10, 15]      # Number of top taxa to display
RANDOM_SEED = 42             # Seed for NMDS reproducibility
```

### Headless Execution

The script uses `matplotlib.use('Agg')` and requires **no display server**. It runs on HPC clusters, Docker containers, and remote servers without modification.

---

## Output Structure

```
results/
├── taxa_barplots/
│   ├── {level}_top{N}_barplot.png           # Grouped barplot (Patient vs Control)
│   ├── {level}_top{N}_per_sample_barplot.png # Per-sample barplot
│   ├── {level}_top{N}_group_means.csv       # Group-level mean abundances
│   └── {level}_top{N}_per_sample_data.csv   # Per-sample abundances (%)
│
├── alpha_diversity/
│   ├── plots/
│   │   └── {level}_alpha_boxplots.png       # Boxplots per diversity metric
│   ├── values/
│   │   └── {level}_alpha_diversity.csv      # Per-sample diversity values
│   └── stats/
│       └── alpha_diversity_stats.csv        # Kruskal–Wallis & Mann–Whitney U
│
└── beta_diversity/
    ├── distance_matrices/
    │   └── {metric}_distance_matrix.csv     # Full pairwise distance matrix
    ├── pcoa/
    │   ├── PCoA_{metric}_top{N}.png         # PCoA scatter with 95% ellipses
    │   └── PCoA_{metric}_coords.csv         # PCoA coordinates
    ├── nmds/
    │   ├── NMDS_{metric}_top{N}.png         # NMDS scatter with 95% ellipses
    │   └── NMDS_{metric}_coords.csv         # NMDS coordinates
    ├── boxplots/
    │   ├── {metric}_boxplot.png             # Within vs between group distances
    │   └── {metric}_boxplot_data.csv
    ├── permanova/
    │   └── PERMANOVA_{metric}.csv           # Per-metric PERMANOVA result
    ├── permdisp/
    │   └── PERMDISP_{metric}.csv            # Per-metric PERMDISP result
    └── beta_stats_summary.csv               # Consolidated statistics
```

Where:
- `{level}` ∈ {`phylum`, `family`, `genus`, `species`}
- `{metric}` ∈ {`Bray_Curtis`, `Jaccard`, `Aitchison`}
- `{N}` ∈ {`10`, `15`}

---

## Methods

### Taxonomic Profiling

Relative abundances are visualized as horizontal stacked barplots. "Unassigned" and "Other" categories are excluded from visualization; the remaining taxa are renormalized to sum to 100%. Taxa are ranked by global mean relative abundance and the top *N* are displayed with a fixed color mapping across all plots.

### Alpha Diversity

Three indices are computed per sample from relative abundance vectors:

- **Shannon index** [[1]](#ref1): $H' = -\sum_{i=1}^{S} p_i \ln(p_i)$
- **Simpson index** [[2]](#ref2): $D = \sum_{i=1}^{S} p_i^2$
- **Inverse Simpson**: $D^{-1} = 1 / D$

Statistical comparison between groups uses the Kruskal–Wallis *H*-test [[3]](#ref3) and pairwise Mann–Whitney *U*-test [[4]](#ref4).

### Beta Diversity

Three distance metrics are computed from the genus-level table:

| Metric | Input | Method |
|--------|-------|--------|
| **Bray–Curtis** [[5]](#ref5) | Relative abundance | $BC_{jk} = 1 - \frac{2 \sum \min(x_{ij}, x_{ik})}{\sum x_{ij} + \sum x_{ik}}$ |
| **Jaccard** [[6]](#ref6) | Presence/absence | $J = 1 - \frac{|A \cap B|}{|A \cup B|}$ |
| **Aitchison** [[7]](#ref7) | CLR-transformed | Euclidean distance on CLR values (pseudocount $10^{-6}$) |

### Ordination

- **PCoA** (Principal Coordinates Analysis) [[8]](#ref8): classical multidimensional scaling on the distance matrix
- **NMDS** (Non-metric Multidimensional Scaling) [[9]](#ref9): via `sklearn.manifold.MDS` with `random_state=42`

95% confidence ellipses are drawn per group.

### Statistical Tests

- **PERMANOVA** [[10]](#ref10): Permutational Multivariate Analysis of Variance with 999 permutations, testing whether group centroids differ in multivariate space
- **PERMDISP** [[11]](#ref11): Permutational test of Multivariate Dispersions, testing homogeneity of within-group dispersions

---

## References

<a id="ref1"></a>[1] Shannon, C. E. (1948). A mathematical theory of communication. *The Bell System Technical Journal*, 27(3), 379–423. https://doi.org/10.1002/j.1538-7305.1948.tb01338.x

<a id="ref2"></a>[2] Simpson, E. H. (1949). Measurement of diversity. *Nature*, 163(4148), 688. https://doi.org/10.1038/163688a0

<a id="ref3"></a>[3] Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion variance analysis. *Journal of the American Statistical Association*, 47(260), 583–621. https://doi.org/10.1080/01621459.1952.10483441

<a id="ref4"></a>[4] Mann, H. B., & Whitney, D. R. (1947). On a test of whether one of two random variables is stochastically larger than the other. *The Annals of Mathematical Statistics*, 18(1), 50–60. https://doi.org/10.1214/aoms/1177730491

<a id="ref5"></a>[5] Bray, J. R., & Curtis, J. T. (1957). An ordination of the upland forest communities of southern Wisconsin. *Ecological Monographs*, 27(4), 325–349. https://doi.org/10.2307/1942268

<a id="ref6"></a>[6] Jaccard, P. (1912). The distribution of the flora in the alpine zone. *New Phytologist*, 11(2), 37–50. https://doi.org/10.1111/j.1469-8137.1912.tb05611.x

<a id="ref7"></a>[7] Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*. Chapman and Hall. https://doi.org/10.1007/978-94-009-4109-0

<a id="ref8"></a>[8] Gower, J. C. (1966). Some distance properties of latent root and vector methods used in multivariate analysis. *Biometrika*, 53(3–4), 325–338. https://doi.org/10.1093/biomet/53.3-4.325

<a id="ref9"></a>[9] Kruskal, J. B. (1964). Multidimensional scaling by optimizing goodness of fit to a nonmetric hypothesis. *Psychometrika*, 29(1), 1–27. https://doi.org/10.1007/BF02289565

<a id="ref10"></a>[10] Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology*, 26(1), 32–46. https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x

<a id="ref11"></a>[11] Anderson, M. J. (2006). Distance-based tests for homogeneity of multivariate dispersions. *Biometrics*, 62(1), 245–253. https://doi.org/10.1111/j.1541-0420.2005.00440.x

### Software

- **NumPy**: Harris, C. R., et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362. https://doi.org/10.1038/s41586-020-2649-2
- **pandas**: McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 56–61. https://doi.org/10.25080/Majora-92bf1922-00a
- **Matplotlib**: Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9(3), 90–95. https://doi.org/10.1109/MCSE.2007.55
- **seaborn**: Waskom, M. L. (2021). seaborn: statistical data visualization. *Journal of Open Source Software*, 6(60), 3021. https://doi.org/10.21105/joss.03021
- **SciPy**: Virtanen, P., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272. https://doi.org/10.1038/s41592-019-0686-2
- **scikit-learn**: Pedregosa, F., et al. (2011). Scikit-learn: Machine Learning in Python. *JMLR*, 12, 2825–2830. https://jmlr.org/papers/v12/pedregosa11a.html

---

## Citation

If you use this pipeline in your research, please cite the relevant methods and software listed in the [References](#references) section above.

---

## License

This project is licensed under the MIT License.
