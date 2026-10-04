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

CONFIG = yaml.safe_load(open("config.yaml"))

import numpy as np
import os
from sklearn.decomposition import PCA


# ---------- BOTH ----------
def load_data(counts_path, metadata_path):
    """Return (counts DataFrame genes x samples, metadata DataFrame indexed by sample).
    Make sure the sample order in counts matches metadata."""

    counts = pd.read_csv(counts_path, index_col=0)
    metadata = pd.read_csv(metadata_path, index_col=0)

    counts = counts.loc[:, metadata.index]

    return counts, metadata

# ---------- Student A ----------
def normalize_log_cpm(counts):
    """log2(counts per million + 1). Return a DataFrame with the same shape."""

    library_sizes = counts.sum(axis=0)
    cpm = counts.div(library_sizes, axis=1) * 1_000_000
    logcpm = np.log2(cpm + 1)

    return logcpm

def plot_pca(logcpm, metadata, out="results/pca.png"):
    """PCA on the samples (transpose!), scatter PC1 vs PC2 coloured by treatment.
    Hint: sklearn.decomposition.PCA on the n_top_features most variable genes."""
    
    n_top = CONFIG.get("n_top_features", 500)

    top_genes = logcpm.var(axis=1).nlargest(n_top).index
    pca_data = logcpm.loc[top_genes].T

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(pca_data)

    os.makedirs(os.path.dirname(out), exist_ok=True)

    for treatment in metadata["treatment"].unique():
        mask = metadata["treatment"] == treatment
        plt.scatter(
            pcs[mask.to_numpy(), 0],
            pcs[mask.to_numpy(), 1],
            label=treatment
        )

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out)
    plt.close()

# ---------- Student B ----------
def run_deseq(counts, metadata):
    """Differential expression treated vs control with pydeseq2.
    Return the results DataFrame (log2FoldChange, padj, ...)."""
    raise NotImplementedError


def top_genes(results, n, out="results/top_genes.csv"):
    """Return and save the n genes with the smallest padj."""
    raise NotImplementedError

# ---------- BOTH ----------
def summary_sentence(results):
    """One sentence answering the research question, e.g. how many genes have padj < 0.05
    and which is the strongest. Return a str."""
    
    significant = results[results["padj"] < 0.05].dropna(subset=["padj"])

    if len(significant) == 0:
        return "No genes were significantly differentially expressed at padj < 0.05."

    strongest_gene = significant["log2FoldChange"].abs().idxmax()

    return (
        f"{len(significant)} genes responded significantly to dexamethasone "
        f"(padj < 0.05); the strongest response was observed for {strongest_gene}."
    )

def main():
    counts, metadata = load_data("data/airway_counts.csv", "data/airway_metadata.csv")
    # Student A: normalize + PCA
    # Student B: DESeq + top genes
    # After the merge: both, then print(summary_sentence(results))
    logcpm = normalize_log_cpm(counts)
    plot_pca(logcpm, metadata)


if __name__ == "__main__":
    main()
