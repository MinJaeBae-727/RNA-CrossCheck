"""Step 5 - volcano and PCA, and a cross-check of the PCA against R.

DESeq2's plotPCA takes the 500 most variable genes of the variance
stabilising transform and runs a principal component analysis on them.
This reproduces that and compares the coordinates with R's.

The sign of a principal component is arbitrary: a PCA that flips an axis
describes the same geometry. The check therefore aligns the sign first.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REF = Path("data/r_reference")
OUT = Path("results")
OUT.mkdir(exist_ok=True)

res = pd.read_csv(OUT / "python_results.csv", index_col="gene_id")
sheet = pd.read_csv(REF / "sample_sheet.csv").set_index("sample")

# ---- volcano -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.4, 5.2), dpi=150)
x = res["log2FoldChange"]
y = -np.log10(res["padj"].clip(lower=1e-300))
sig = (res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)

ax.scatter(x[~sig], y[~sig], s=5, alpha=0.3, color="#b0b0b0", edgecolors="none",
           label=f"not significant ({(~sig).sum()})")
ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=6, alpha=0.7, color="#d62728",
           edgecolors="none", label=f"up ({(sig & (x > 0)).sum()})")
ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=6, alpha=0.7, color="#1f77b4",
           edgecolors="none", label=f"down ({(sig & (x < 0)).sum()})")
ax.axhline(-np.log10(0.05), color="#444", lw=0.6, ls="--")
ax.axvline(1, color="#444", lw=0.6, ls="--")
ax.axvline(-1, color="#444", lw=0.6, ls="--")
ax.set_xlabel("log2 fold change")
ax.set_ylabel("-log10 adjusted P")
ax.set_title("WT_Cold vs WT_NoCold (PyDESeq2)", fontsize=11)
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "volcano.png", bbox_inches="tight")
print(f"volcano written to {OUT / 'volcano.png'}")

# ---- PCA on the 500 most variable VST genes ------------------------------
# Read the VST matrix R produced, so the comparison isolates the PCA itself.
vst = pd.read_csv(REF / "vst.csv", index_col="gene_id")

# R: rv <- rowVars(assay(vsd)); select <- order(rv, decreasing = TRUE)[1:500]
rv = vst.var(axis=1, ddof=1)
top = rv.sort_values(ascending=False).index[:500]
mat = vst.loc[top].T                       # samples in rows

# R: prcomp(t(assay(vsd)[select, ]))  - centres columns, does not scale
centred = mat - mat.mean(axis=0)
u, s, vt = np.linalg.svd(centred, full_matrices=False)
scores = u * s                              # principal component scores
var_explained = s**2 / np.sum(s**2)

py_pca = pd.DataFrame({"PC1": scores[:, 0], "PC2": scores[:, 1]}, index=mat.index)

r_pca = pd.read_csv(REF / "pca.csv", index_col="sample")
py_pca = py_pca.loc[r_pca.index]

print("\nPCA cross-check")
for pc in ["PC1", "PC2"]:
    a, b = r_pca[pc], py_pca[pc]
    # Align the arbitrary sign before comparing.
    if np.corrcoef(a, b)[0, 1] < 0:
        b = -b
        py_pca[pc] = b
    print(f"  {pc}  Pearson r {stats.pearsonr(a, b)[0]:.10f}"
          f"   max |difference| {(a - b).abs().max():.3e}")

print(f"  variance explained, Python : PC1 {var_explained[0]*100:.0f}%  PC2 {var_explained[1]*100:.0f}%")
print(f"  variance explained, R      : PC1 {r_pca['pc1_percent'].iloc[0]:.0f}%  PC2 {r_pca['pc2_percent'].iloc[0]:.0f}%")

fig, ax = plt.subplots(figsize=(5.6, 4.8), dpi=150)
colors = {"WT_NoCold": "#1f77b4", "WT_Cold": "#d62728",
          "mutant_NoCold": "#7fb3d5", "mutant_Cold": "#e8928c"}
for g, sub in py_pca.groupby(sheet.loc[py_pca.index, "group"]):
    ax.scatter(sub["PC1"], sub["PC2"], s=55, label=g, color=colors.get(g), edgecolors="white", lw=0.8)
ax.set_xlabel(f"PC1  ({var_explained[0]*100:.0f} % variance)")
ax.set_ylabel(f"PC2  ({var_explained[1]*100:.0f} % variance)")
ax.set_title("PCA of the 500 most variable genes (VST)", fontsize=11)
ax.legend(frameon=False, fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "pca.png", bbox_inches="tight")
print(f"\nPCA figure written to {OUT / 'pca.png'}")
