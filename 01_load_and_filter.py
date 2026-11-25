"""Step 1 - load the counts and reproduce the R pre-filter.

Target: DESeq2 kept 5857 of 6000 genes. This script must print the same number.
Each comment starting with "R:" shows the line this replaces.
"""

from pathlib import Path
import pandas as pd

REF = Path("data/r_reference")

# R: counts <- read.csv("counts.csv", row.names = 1, check.names = FALSE)
counts = pd.read_csv(REF / "counts.csv", index_col="gene_id")

# R: sheet <- read.csv("sample_sheet.csv")
sheet = pd.read_csv(REF / "sample_sheet.csv")

print("counts:", counts.shape)        # R: dim(counts)
print(counts.iloc[:3, :4])            # R: head(counts[, 1:4], 3)
print()
print(sheet.head(4))                  # R: head(sheet, 4)

# The filter inside run_deseq():
#   min_samples <- min(table(sheet$group))
#   keep        <- rowSums(counts >= 10) >= min_samples
min_count = 10
min_samples = sheet["group"].value_counts().min()

# (counts >= 10) is a True/False table; .sum(axis=1) counts True per row.
# R: rowSums(counts >= min_count)
n_passing = (counts >= min_count).sum(axis=1)
keep = n_passing >= min_samples

filtered = counts[keep]               # R: counts[keep, ]

print()
print(f"min_count   = {min_count}")
print(f"min_samples = {min_samples}")
print(f"kept {keep.sum()} of {len(counts)} genes")
print(f"filtered: {filtered.shape}")
