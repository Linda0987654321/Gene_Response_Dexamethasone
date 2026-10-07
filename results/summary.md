# Transcriptomic analysis summary

## Research question
Which genes respond to dexamethasone in airway smooth muscle cells?

## Overview
This analysis used `analysis.py` to compare gene expression between 4 control and 4 dexamethasone-treated airway smooth muscle cell samples, using raw counts for 19,271 genes. For the PCA, counts were converted to counts per million (CPM) for each sample and log2-transformed as `log2(CPM + 1)`. The script produces a PCA plot.

In addition, it creates a table of the 20 genes with the smallest adjusted p-values and prints a summary of significant genes.

## Results and Interpretation
We found 2,201 genes with significantly different expression between dexamethasone-treated and control cells (adjusted p-value < 0.05). Some showed higher expression after treatment, while others showed lower expression.
ENSG00000152583 had the smallest adjusted p-value (4.48 × 10⁻⁷¹), providing the strongest statistical evidence for a difference, although not necessarily the largest change in expression.

Treated samples (orange) appear below control samples (blue) at similar positions along PC1. This suggests a treatment-related difference in gene expression along PC2. However, the groups do not form two clearly separate clusters, and both vary along PC1, suggesting other differences between samples. Lower PC2 values do not mean lower overall gene expression.

## Output files
- `results/pca.png`: a plot comparing the samples.
- `results/top_genes.csv`: a table of the top-ranked genes, with expression changes and adjusted p-values.

The PCA plot provides an overview of the samples, while the differential expression analysis identifies individual genes associated with the treatment.
