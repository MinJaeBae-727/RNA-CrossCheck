# rnaseq-crosscheck

The same RNA-seq differential expression analysis, implemented twice -
once in R (DESeq2) and once in Python (PyDESeq2) - and checked against
each other gene by gene.

## Why

Two independent implementations of one analysis should agree. Where they
do not, one of them is wrong. This repository makes that check explicit
and reproducible.

## Data

Simulated counts from `make_demo_rnaseq()` in PhytoOmics Studio:
6,000 Arabidopsis genes, 2 genotypes x 2 treatments x 3 replicates,
negative binomial with a known ground truth. No measured data is used,
so the repository can stay public.

## Layout

    export_reference.R        R side: simulate, run DESeq2, write CSVs
    data/r_reference/         the R results Python has to reproduce
    01_load_and_filter.py     load counts, reproduce the pre-filter

## Status

- [x] Step 1 - pre-filter: 5,857 of 6,000 genes kept, same as R
- [ ] Step 2 - size factors
- [ ] Step 3 - differential expression with PyDESeq2
- [ ] Step 4 - compare log2 fold changes and adjusted P values
- [ ] Step 5 - volcano and PCA plots

## Running it

    Rscript export_reference.R
    .venv/Scripts/python 01_load_and_filter.py
