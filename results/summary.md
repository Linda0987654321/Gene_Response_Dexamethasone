# Transcriptomic analysis summary

## Research question
Which genes respond to dexamethasone in airway smooth muscle cells?

## What we did
We analysed gene expression data from eight samples, comparing dexamethasone-treated cells with control cells. We first checked that the counts and sample information matched.

Student A converted the counts to counts per million (CPM), applied a log transformation, and created a principal component analysis (PCA) plot. This plot helps us explore similarities and differences between samples based on their gene expression. Each point represents one sample, coloured by treatment. Samples that are close together have more similar expression patterns, while samples farther apart have more different patterns.

Student B used PyDESeq2 with the original counts to identify genes with different expression between treated and control cells. The results show the size and direction of each change, together with an adjusted p-value. Genes were ranked by adjusted p-value.

## Results
We found 2,201 genes with significantly different expression between dexamethasone-treated and control cells (adjusted p-value < 0.05). Some showed higher expression after treatment, while others showed lower expression.
ENSG00000152583 had the smallest adjusted p-value (4.48 × 10⁻⁷¹), providing the strongest statistical evidence for a difference, although not necessarily the largest change in expression.

Treated samples (orange) appear below control samples (blue) at similar positions along PC1. This suggests a treatment-related difference in gene expression along PC2. However, the groups do not form two clearly separate clusters, and both vary along PC1, suggesting other differences between samples. Lower PC2 values do not mean lower overall gene expression.

## Output files
- `results/pca.png`: a plot comparing the samples.
- `results/top_genes.csv`: a table of the top-ranked genes, with expression changes and adjusted p-values.

The PCA plot provides an overview of the samples, while the differential expression analysis identifies individual genes associated with the treatment.
