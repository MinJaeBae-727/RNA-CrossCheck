# Export the R-side reference: simulated counts, sample sheet and DESeq2 results.
# Python must reproduce the numbers in data/r_reference/deseq2_results.csv.
suppressPackageStartupMessages(pkgload::load_all("E:/DATA/PhytoOmicsStudio", quiet = TRUE))
set.seed(1)
demo <- make_demo_rnaseq(n_genes = 6000, seed = 1, species = "ath")

out <- "data/r_reference"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

counts <- demo$counts
write.csv(data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
          file.path(out, "counts.csv"), row.names = FALSE)
write.csv(demo$sheet, file.path(out, "sample_sheet.csv"), row.names = FALSE)

dds <- run_deseq(demo, demo$sheet, factor = "group", min_count = 10)
res <- get_de_results(dds, numerator = "WT_Cold", denominator = "WT_NoCold",
                      alpha = 0.05, shrink = "none")
write.csv(res, file.path(out, "deseq2_results.csv"), row.names = FALSE)

sf <- DESeq2::sizeFactors(dds)
write.csv(data.frame(sample = names(sf), size_factor = as.numeric(sf)),
          file.path(out, "size_factors.csv"), row.names = FALSE)

writeLines(c(
  sprintf("R              %s", getRversion()),
  sprintf("DESeq2         %s", packageVersion("DESeq2")),
  sprintf("genes_before   %d", dds@metadata$filter$n_before),
  sprintf("genes_after    %d", dds@metadata$filter$n_after),
  sprintf("design         %s", dds@metadata$filter$design),
  sprintf("contrast       WT_Cold vs WT_NoCold")
), file.path(out, "run_info.txt"))

cat("counts     ", nrow(counts), "x", ncol(counts), "\n")
cat("after filter", dds@metadata$filter$n_after, "genes\n")
cat("results    ", nrow(res), "rows,", ncol(res), "cols\n")
cat("columns    ", paste(names(res), collapse = ", "), "\n")
