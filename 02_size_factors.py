"""Step 2 - compute DESeq2 size factors (median of ratios) in Python.

DESeq2 normalises for sequencing depth with one factor per sample:

    loggeomeans_i = mean over samples of log(count_ij)        (per gene i)
    sizefactor_j  = exp( median_i [ log(count_ij) - loggeomeans_i ] )

The median runs over genes whose log geometric mean is finite (no zero in
any sample) and whose count in sample j is above zero.

Target: data/r_reference/size_factors.csv, written by DESeq2.
"""

from pathlib import Path
import numpy as np
import pandas as pd

REF = Path("data/r_reference")

counts = pd.read_csv(REF / "counts.csv", index_col="gene_id")
sheet = pd.read_csv(REF / "sample_sheet.csv")

# Same pre-filter as step 1: DESeq2 computes size factors after filtering.
min_samples = sheet["group"].value_counts().min()
counts = counts[(counts >= 10).sum(axis=1) >= min_samples]

# np.log(0) is -inf and raises a warning; errstate silences it because the
# -inf values are exactly what the finite-check below removes.
# R: log_counts <- log(counts)
with np.errstate(divide="ignore"):
    log_counts = np.log(counts)

# R: loggeomeans <- rowMeans(log(counts))
loggeomeans = log_counts.mean(axis=1)

# R: finite <- is.finite(loggeomeans)
finite = np.isfinite(loggeomeans)

size_factors = {}
for sample in counts.columns:
    # R: use <- finite & counts[, j] > 0
    use = finite & (counts[sample] > 0)
    # R: exp(median((log(counts[, j]) - loggeomeans)[use]))
    size_factors[sample] = np.exp(np.median(log_counts.loc[use, sample] - loggeomeans[use]))

python_sf = pd.Series(size_factors, name="python")

# ---- compare with R -------------------------------------------------------
r_sf = pd.read_csv(REF / "size_factors.csv", index_col="sample")["size_factor"]
r_sf.name = "R"

cmp = pd.concat([r_sf, python_sf], axis=1)
cmp["difference"] = cmp["python"] - cmp["R"]
cmp["rel_error"] = cmp["difference"] / cmp["R"]

pd.set_option("display.float_format", lambda v: f"{v: .10f}")
print(cmp)
print()
print(f"largest absolute difference : {cmp['difference'].abs().max():.3e}")
print(f"largest relative error      : {cmp['rel_error'].abs().max():.3e}")

# np.allclose asks "are these the same to within rounding?"
if np.allclose(cmp["python"], cmp["R"], rtol=1e-10, atol=0):
    print("\nMATCH - Python and R size factors agree to 1e-10 relative.")
else:
    print("\nMISMATCH - the two implementations disagree.")
