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


# ---------- BOTH ----------
def load_data(counts_path, metadata_path):
    """Return (counts DataFrame genes x samples, metadata DataFrame indexed by sample).
    Make sure the sample order in counts matches metadata."""
    raise NotImplementedError


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
    raise NotImplementedError


def top_genes(results, n, out="results/top_genes.csv"):
    """Return and save the n genes with the smallest padj."""
    raise NotImplementedError


# ---------- BOTH ----------
def summary_sentence(results):
    """One sentence answering the research question, e.g. how many genes have padj < 0.05
    and which is the strongest. Return a str."""
    raise NotImplementedError


def main():
    counts, metadata = load_data("data/airway_counts.csv", "data/airway_metadata.csv")
    # Student A: normalize + PCA
    # Student B: DESeq + top genes
    # After the merge: both, then print(summary_sentence(results))
    raise NotImplementedError


if __name__ == "__main__":
    main()
