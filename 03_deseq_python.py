"""Step 3 - run the differential expression analysis with PyDESeq2.

Mirrors run_deseq() + get_de_results() from PhytoOmics Studio:
    design    ~ group
    filter    at least 10 reads in at least 3 samples
    contrast  WT_Cold vs WT_NoCold

Writes results/python_results.csv.
"""

from pathlib import Path
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

REF = Path("data/r_reference")
OUT = Path("results")
OUT.mkdir(exist_ok=True)

counts = pd.read_csv(REF / "counts.csv", index_col="gene_id")
sheet = pd.read_csv(REF / "sample_sheet.csv").set_index("sample")

# Same pre-filter as R.
min_samples = sheet["group"].value_counts().min()
counts = counts[(counts >= 10).sum(axis=1) >= min_samples]
print(f"{counts.shape[0]} genes after filtering, {counts.shape[1]} samples")

# PyDESeq2 wants samples in rows and genes in columns - the transpose of the
# R convention. R: t(counts)
counts_t = counts.T
# Keep the sample order of the count table.
metadata = sheet.loc[counts_t.index]

# R: DESeqDataSetFromMatrix(counts, cd, ~ group); DESeq(dds)
dds = DeseqDataSet(counts=counts_t, metadata=metadata, design="~group", quiet=True)
dds.deseq2()

# R: results(dds, contrast = c("group", "WT_Cold", "WT_NoCold"), alpha = 0.05)
stat = DeseqStats(dds, contrast=["group", "WT_Cold", "WT_NoCold"], alpha=0.05, quiet=True)
stat.summary()
res = stat.results_df

res.index.name = "gene_id"
res = res.sort_values("pvalue")
res.to_csv(OUT / "python_results.csv")

print()
print(res.head(8))
print()
print(f"tested          {len(res)} genes")
print(f"padj < 0.05     {(res['padj'] < 0.05).sum()}")
print(f"padj is NA      {res['padj'].isna().sum()}")
print(f"written to      {OUT / 'python_results.csv'}")
