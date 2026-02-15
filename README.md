# Interpretable Rule-Based Prediction of Cisplatin Response in Lung Cancer

## Overview
This project applies **R.ROSETTA**, an interpretable rule-based machine learning framework based on rough set theory, to predict cisplatin drug response in lung cancer cell lines using gene expression data from the GDSC2 database.

Unlike black-box models (neural networks, SVM, random forests), R.ROSETTA produces **human-readable IF-THEN rules** that reveal co-predictive gene mechanisms underlying drug sensitivity and resistance.

## Key Results

| Metric | Value |
|--------|-------|
| Accuracy (10-fold CV) | **73.9%** |
| Permutation p-value | **0** (n=100) |
| Total rules discovered | **241** |
| Samples | 127 lung cancer cell lines |
| Features | 43 genes (LASSO + Boruta selection) |
| Classes | Sensitive (n=64) vs Resistant (n=63) |

## Key Findings

### Biologically Validated Discoveries
- **EGFR** (21 resistance rules): Known oncogene; EGFR-high lung cancers show platinum resistance
- **SLFN11** (sensitivity rules): Established cisplatin sensitivity biomarker in literature
- **MSH2** (sensitivity rules): DNA mismatch repair gene; functional MMR enhances cisplatin efficacy

### Novel Findings
- **TUFT1** is the top hub gene (35 rules) connecting both resistance and sensitivity networks — a potential novel cisplatin response biomarker
- **SERINC3** (30 rules) and **PRPF40A** (20 rules) emerge as key co-predictive features not previously linked to cisplatin response

### Example Rules
```
IF WASF3=high AND PRPF40A=high THEN Sensitive (accuracy=100%, support=20)
IF SSBP1=low AND PRPF40A=low THEN Resistant (accuracy=100%, support=19)
IF EGFR=high AND STOML2=low THEN Resistant (accuracy=100%, support=16)
IF SLFN11=high AND PRPF40A=high THEN Sensitive (accuracy=100%, support=13)
```

## Pipeline (Following Komorowski RBM Guidelines)

1. **Biological Question**: Predict cisplatin sensitivity in lung cancer
2. **Data**: GDSC2 via PharmacoGx R package
3. **Preprocessing**: NA imputation, variance filtering, z-score normalization
4. **Feature Selection**: LASSO (34 genes) + Boruta (11 genes) → 43 gene union
5. **Discretization**: Equal frequency binning (low/medium/high)
6. **RBM**: R.ROSETTA with 10-fold CV, StandardVoter classifier
7. **Validation**: Permutation test (p=0)
8. **Post-processing**: Rule network analysis, hub gene identification

## Files

| File | Description |
|------|-------------|
| `analysis.R` | Complete reproducible R script |
| `cisplatin_rules.csv` | All 241 rules with statistics |
| `expression_data.csv` | Normalized gene expression (127 x 43) |
| `response_labels.csv` | Sensitive/Resistant labels |
| `lasso_features.csv` | LASSO-selected genes |
| `boruta_features.csv` | Boruta-selected genes |
| `network_edges.csv` | Rule co-occurrence network edges |
| `figures.pdf` | PCA, distribution, permutation test, network plots |
| `MODEL_SUMMARY.txt` | Full model summary |

## Requirements
```r
install.packages(c("BiocManager", "glmnet", "Boruta", "igraph", "devtools"))
BiocManager::install("PharmacoGx")
devtools::install_github("komorowskilab/R.ROSETTA")
```

## How to Reproduce
```r
source("analysis.R")
```
Note: Initial data download (~1.5 GB) takes 15-30 minutes.

## References

- Komorowski, J. "Developing RBMs Guidelines" - Knowledge-based Systems in Bioinformatics
- Dramiński, M. et al. "R.ROSETTA: an interpretable machine learning framework" (2021) BMC Bioinformatics
- Yang, W. et al. "Genomics of Drug Sensitivity in Cancer (GDSC)" Nucleic Acids Research (2013)

## License
MIT License
