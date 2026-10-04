# Group 1 - Transcriptomics: glucocorticoid response in airway smooth muscle

**Research question:** Which genes respond to dexamethasone in airway smooth muscle cells?

**Data:** `data/airway_counts.csv` (19,271 genes x 8 samples, raw counts, GEO GSE52778) and `data/airway_metadata.csv` (treatment, cell line).

Everything happens in **one file, `analysis.py`**. Each function is a stub that raises
`NotImplementedError`. The file header says who implements what.

| Owner | Functions |
|---|---|
| Student A | `normalize_log_cpm`, `plot_pca` |
| Student B | `run_deseq`, `top_genes` |
| **Both** (this is where the merge conflict happens) | `load_data`, `main`, `summary_sentence` |

Shared parameters live in `config.yaml`. Both students must set the value marked
`# BOTH` — you will disagree, and Git cannot decide for you.

Run with `python analysis.py`. It must run without error after your merge.

## Results (fill in after the merge)

<!-- one sentence answering the research question, one figure -->

## Reflection

1. What caused each merge conflict?
2. How could branching strategy or file layout have avoided it?
3. What is the difference between the history produced by `git pull` and `git pull --rebase`?
