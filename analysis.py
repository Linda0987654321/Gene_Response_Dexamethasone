"""
Shared analysis script. Both students edit THIS file on their own branches.

Student A implements : normalize_log_cpm, plot_pca
Student B implements : run_deseq, top_genes
BOTH implement       : load_data, main, summary_sentence   <- expect a merge conflict here

Run: python analysis.py
"""
import yaml
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import math
import os

CONFIG = yaml.safe_load(open("config.yaml"))

from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ---------- BOTH ----------
def load_data(counts_path, metadata_path):
    """Return (counts DataFrame genes x samples, metadata DataFrame indexed by sample).
    Make sure the sample order in counts matches metadata."""

    # The first column contains gene IDs, so use it as the row index.
    counts = pd.read_csv(counts_path, index_col=0)

    # Use the sample column as the metadata row index.
    # This lets us look up metadata by sample ID.
    metadata = pd.read_csv(metadata_path, index_col="sample")

    # Keep samples that appear in both files, preserving the order in counts.
    #samples = [sample for sample in counts.columns if sample in metadata.index]

    #if not samples:
        #raise ValueError("No sample IDs match between counts and metadata.")

    # Check that every sample in the counts file has metadata.
    missing_metadata = [
        sample for sample in counts.columns
        if sample not in metadata.index
    ]

    # Check that every metadata sample has counts.
    missing_counts = [
        sample for sample in metadata.index
        if sample not in counts.columns
    ]

    if missing_metadata or missing_counts:
        raise ValueError(
            f"Sample IDs do not match. "
            f"Missing from metadata: {missing_metadata}; "
            f"missing from counts: {missing_counts}"
        )

    # All IDs match. Reorder metadata to match the counts columns.
    metadata = metadata.loc[counts.columns]
       
    # Raw read counts cannot be negative.
    if (counts < 0).any().any():
        raise ValueError("Counts must be non-negative.")

    # PyDESeq2 requires integer raw counts.
    if not (counts.to_numpy() == counts.to_numpy().astype(int)).all():
        raise ValueError("Counts must be integers for differential expression.")

    return counts.astype(int), metadata


# ---------- Student A ----------
def normalize_log_cpm(counts):
    """log2(counts per million + 1). Return a DataFrame with the same shape."""
    raise NotImplementedError


def plot_pca(logcpm, metadata, out="results/pca.png"):
    """PCA on the samples (transpose!), scatter PC1 vs PC2 coloured by treatment.
    Hint: sklearn.decomposition.PCA on the n_top_features most variable genes."""
    raise NotImplementedError


# ---------- Student B ----------
def run_deseq(counts, metadata):
    """Differential expression treated vs control with pydeseq2.
    Return the results DataFrame (log2FoldChange, padj, ...)."""

    # Test each gene for a difference between treated and control samples.

    # PyDESeq2 expects samples in rows and genes in columns.
    # Transposing the table changes the index from gene ID to sample ID
    count_data = counts.T

    # Select rows from metadata and match them with the sample IDs 
    # in the transposed counts table, keeping only the "treatment" column
    sample_metadata = metadata.loc[count_data.index, ["treatment"]]

    # Creating a DESeq2 dataset object from the counts and sample information, 
    # to be used for differential expression analysis.
    dds = DeseqDataSet(
        counts=count_data,
        metadata=sample_metadata,   
        design="~treatment", #gene counts should be analyzed in relation to the treatment
    )

    # Run the differential expression analysis.
    # PyDESeq2 performs the statistical comparison of gene counts, estimating 
    # both the size and direction of each difference and the statistical evidence for it.
    dds.deseq2()

    # Compare dexamethasone samples against control samples.
    # Setting up statistical comparison using the fitted dds dataset.
    # Positive log2FoldChange means higher expression under dexamethasone.
    # Negative log2FoldChange means lower expression under dexamethasone.
    stats = DeseqStats(
        dds,
        contrast=["treatment", "treated", "control"],
    )
    stats.summary()

    # Copy the results table, which includes fold changes and adjusted p-values.
    results = stats.results_df.copy()
    results.index.name = "gene_id"
    return results


def top_genes(results, n, out="results/top_genes.csv"):
    """Return and save the n genes with the smallest padj."""

    #Save a ranked selection of the genes found with PyDESeq
    #padj is the adjusted p-value provided in the result table from run_deseq

    if "padj" not in results.columns:
        raise ValueError("Results must contain a 'padj' column.")

    n = int(n)
    if n <= 0:
        raise ValueError("n must be positive.")

    # Ignore genes without an adjusted p-value.
    # Sort by adjusted p-value, then keep the first n genes.
    top = results.dropna(subset=["padj"]).sort_values("padj").head(n)

    # Create the results folder if it does not already exist.
    os.makedirs(os.path.dirname(out), exist_ok=True)

    # Save the table with gene IDs as the first column.
    top.to_csv(out, index=True, index_label="gene_id")
    return top

# ---------- BOTH ----------
def summary_sentence(results):
    """One sentence answering the research question, e.g. how many genes have padj < 0.05
    and which is the strongest. Return a str."""

    # Read the adjusted p-value cutoff from config.yaml.
    threshold = float(CONFIG["padj_threshold"])

    # Ignore genes without an adjusted p-value.
    tested = results.dropna(subset=["padj"])

    if tested.empty:
        return "No genes had an adjusted p-value available for testing."

    # Count genes whose adjusted p-value is below the cutoff.
    significant = tested[tested["padj"] < threshold]

    # Find the gene with the smallest adjusted p-value.
    strongest_gene = tested["padj"].idxmin()
    strongest_padj = tested.loc[strongest_gene, "padj"]

    # Return one sentence summarizing the results.
    return (
        f"{len(significant)} genes had padj < {threshold:g}; "
        f"the strongest result was {strongest_gene} "
        f"(padj={strongest_padj:.3g})."
    )


def main():
    # Load the file paths from config.yaml and read both data tables.
    counts, metadata = load_data(
        CONFIG["counts"],
        CONFIG["metadata"],
    )

    # Student A: normalize the counts and create the PCA plot.
    logcpm = normalize_log_cpm(counts)
    plot_pca(logcpm, metadata)

    # Student B: run differential expression and save the top genes.
    results = run_deseq(counts, metadata)
    top_genes(results, CONFIG["n_top_features"])

    # BOTH: print the summary sentence.
    print(summary_sentence(results))


if __name__ == "__main__":
    # Run main() only when this file is executed directly,
    # for example with: python analysis.py
    main()
