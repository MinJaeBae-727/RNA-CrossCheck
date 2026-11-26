"""Step 4 - compare the Python results against the R results, gene by gene.

Two independent implementations of one analysis should agree. This script
quantifies how far they do.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")          # write image files, do not open a window
import matplotlib.pyplot as plt

REF = Path("data/r_reference")
OUT = Path("results")
OUT.mkdir(exist_ok=True)

r = pd.read_csv(REF / "deseq2_results.csv", index_col="gene_id")
p = pd.read_csv(OUT / "python_results.csv", index_col="gene_id")

print(f"R      {len(r)} genes")
print(f"Python {len(p)} genes")
print(f"shared {len(r.index.intersection(p.index))} genes")

# R: merge(r, p, by = "gene_id")
m = r.join(p, how="inner", lsuffix="_R", rsuffix="_py") if False else r.add_suffix("_R").join(p.add_suffix("_py"), how="inner")

def agreement(col, label):
    a, b = m[f"{col}_R"], m[f"{col}_py"]
    ok = a.notna() & b.notna()
    a, b = a[ok], b[ok]
    pear = stats.pearsonr(a, b)[0]
    spear = stats.spearmanr(a, b)[0]
    print(f"\n{label}  (n = {ok.sum()})")
    print(f"  Pearson r            {pear:.10f}")
    print(f"  Spearman rho         {spear:.10f}")
    print(f"  max |difference|     {(a - b).abs().max():.3e}")
    print(f"  median |difference|  {(a - b).abs().median():.3e}")
    return pear

print("\n" + "=" * 60)
r_lfc = agreement("log2FoldChange", "log2 fold change")
r_bm = agreement("baseMean", "base mean")

# padj spans many orders of magnitude, so compare it on the log scale.
ok = m["padj_R"].notna() & m["padj_py"].notna()
la, lb = -np.log10(m.loc[ok, "padj_R"] + 1e-300), -np.log10(m.loc[ok, "padj_py"] + 1e-300)
print(f"\n-log10 adjusted P  (n = {ok.sum()})")
print(f"  Pearson r            {stats.pearsonr(la, lb)[0]:.10f}")
print(f"  max |difference|     {(la - lb).abs().max():.3e}")

# ---- do the two agree on which genes are significant? --------------------
sig_r = (m["padj_R"] < 0.05).fillna(False)
sig_p = (m["padj_py"] < 0.05).fillna(False)

both = (sig_r & sig_p).sum()
only_r = (sig_r & ~sig_p).sum()
only_p = (~sig_r & sig_p).sum()
neither = (~sig_r & ~sig_p).sum()

print("\n" + "=" * 60)
print("Significant at padj < 0.05")
print(f"  R only      {only_r}")
print(f"  Python only {only_p}")
print(f"  both        {both}")
print(f"  neither     {neither}")
jaccard = both / (both + only_r + only_p) if (both + only_r + only_p) else float("nan")
print(f"  agreement   {(both + neither) / len(m) * 100:.4f} %   Jaccard {jaccard:.6f}")

if only_r or only_p:
    print("\nGenes called by one implementation only:")
    disagree = m[sig_r ^ sig_p][["baseMean_R", "log2FoldChange_R", "padj_R", "padj_py"]]
    print(disagree.sort_values("padj_R").head(10))

# ---- figure --------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(10, 4.6), dpi=150)

ax[0].scatter(m["log2FoldChange_R"], m["log2FoldChange_py"], s=4, alpha=0.3,
              edgecolors="none", color="#1f77b4")
lim = [m[["log2FoldChange_R", "log2FoldChange_py"]].min().min(),
       m[["log2FoldChange_R", "log2FoldChange_py"]].max().max()]
ax[0].plot(lim, lim, color="#d62728", lw=0.8, ls="--")
ax[0].set_xlabel("log2 fold change, R (DESeq2)")
ax[0].set_ylabel("log2 fold change, Python (PyDESeq2)")
ax[0].set_title(f"Pearson r = {r_lfc:.8f}", fontsize=10)

ax[1].scatter(la, lb, s=4, alpha=0.3, edgecolors="none", color="#2ca02c")
lim2 = [0, max(la.max(), lb.max())]
ax[1].plot(lim2, lim2, color="#d62728", lw=0.8, ls="--")
ax[1].set_xlabel("-log10 adjusted P, R")
ax[1].set_ylabel("-log10 adjusted P, Python")
ax[1].set_title("Adjusted P values", fontsize=10)

for a in ax:
    a.spines[["top", "right"]].set_visible(False)
fig.suptitle("R (DESeq2) vs Python (PyDESeq2), WT_Cold vs WT_NoCold", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "r_vs_python.png", bbox_inches="tight")
print(f"\nfigure written to {OUT / 'r_vs_python.png'}")
