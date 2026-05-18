# RNA-seq DEG Analysis: DESeq2 vs Statistical Methods Comparison

A bioinformatics project comparing two approaches to identifying **Differentially Expressed Genes (DEGs)** from RNA-seq data:

- **Python**: Log2-CPM normalization + Student's t-test (implemented in a Jupyter Notebook)
- **R**: DESeq2 (negative binomial model, standard for count-based RNA-seq data)

The dataset comes from a study of **diabetic patients vs. healthy controls**, with samples collected across infection and non-infection conditions.

---

## Project Structure

```
rna-seq-deg-analysis/
├── data/
│   ├── Diabetes_infection_count.csv   # Raw gene count matrix (genes × samples)
│   └── Diabetes_infection_pheno.csv   # Sample metadata (groups, clinical labels)
├── notebooks/
│   └── RNA-seq_DEG_Analysis.ipynb     # Python pipeline: normalization, t-test, heatmap
├── scripts/
│   └── DESeq2_analysis.R              # R pipeline: DESeq2, MA plot, Venn diagram
├── results/                           # Generated outputs (CSVs, plots)
└── README.md
```

---

## Dataset

| File | Description |
|------|-------------|
| `Diabetes_infection_count.csv` | Raw read counts matrix — rows are genes, columns are sample IDs |
| `Diabetes_infection_pheno.csv` | Sample metadata including group labels (`CL4`): `Healthy` or `Infected` |

> **Note:** Sample groups are derived from the `CL4` column in the phenotype file. `Infection` labels are recoded as `Infected` for consistency across both pipelines.

---

## Methods

### Python Pipeline (`notebooks/RNA-seq_DEG_Analysis.ipynb`)

1. **Data loading & cleaning** — fill NAs, convert to integers
2. **Filtering** — remove genes with mean count ≤ 10
3. **Normalization** — Log2 CPM (Counts Per Million) with pseudocount of 1
4. **Group splitting** — separate Healthy vs. Infected samples using metadata
5. **Statistical testing** — independent samples t-test per gene (`scipy.stats.ttest_ind`)
6. **DEG filtering** — p-value < 0.05 and |log2 fold change| > 0.4
7. **Visualization** — expression distribution histogram, clustered heatmap of top 20 DEGs

### R Pipeline (`scripts/DESeq2_analysis.R`)

1. **Data loading** — count matrix + phenotype table
2. **Metadata preparation** — align sample IDs, set `Healthy` as reference
3. **Pre-filtering** — keep genes with ≥ 10 total counts
4. **DESeq2 analysis** — negative binomial model, Wald test
5. **Results** — adjusted p-value (Benjamini-Hochberg), log2 fold change
6. **Visualization** — MA plot, Venn diagram comparing top 20 DEGs from both methods

---

## Comparison

Both pipelines extract the **top 20 DEGs** ranked by p-value, then compare the overlap using a Venn diagram. This highlights the differences in sensitivity and statistical modeling between:

- A simple parametric test on normalized data (t-test)
- A count-aware model designed for RNA-seq (DESeq2)

---

## Requirements

### Python

```
pandas
numpy
scipy
matplotlib
seaborn
```

Install with:
```bash
pip install pandas numpy scipy matplotlib seaborn
```

### R

```r
BiocManager::install("DESeq2")
install.packages(c("tidyverse", "ggVennDiagram"))
```

---

## How to Run

### Python Notebook

```bash
jupyter notebook notebooks/RNA-seq_DEG_Analysis.ipynb
```

### R Script

```bash
Rscript scripts/DESeq2_analysis.R
```

> Run the Python notebook first — it exports `results/DEG_20_ttst.csv` which the R script uses for the Venn diagram comparison.

---

## Results

After running both pipelines, the `results/` folder will contain:

| File | Description |
|------|-------------|
| `top20_DEGs_DESeq2.csv` | Top 20 DEGs from DESeq2 |
| `DEG_20_ttst.csv` | Top 20 DEGs from Python t-test |
| `MA_plot.png` | MA plot from DESeq2 results |
| `venn_diagram.png` | Overlap between DESeq2 and t-test top 20 genes |

---

## Author

**Alfek** — Bioinformatics project, 2025  
Data: Diabetic infection cohort RNA-seq study
