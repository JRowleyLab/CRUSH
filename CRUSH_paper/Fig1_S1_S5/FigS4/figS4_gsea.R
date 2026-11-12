#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(fgsea); library(data.table)})

a <- commandArgs(trailingOnly = TRUE)
r <- fread(a[1], header = FALSE)
stats <- setNames(r$V2, r$V1); stats <- stats[is.finite(stats)]
gmt <- gmtPathways(a[2])
gmt <- lapply(gmt, function(g) unique(toupper(sub(",.*", "", trimws(g)))))
gmt <- lapply(gmt, function(g) intersect(g, names(stats)))
gmt <- gmt[lengths(gmt) >= 5 & lengths(gmt) <= 500]
set.seed(1)
res <- fgsea(gmt, stats, minSize = 5, maxSize = 500, nPermSimple = 10000)
res <- res[order(-abs(NES))]
out <- copy(res); out[, leadingEdge := sapply(leadingEdge, paste, collapse = ",")]
fwrite(out, a[3], sep = "\t")
cat(sprintf("%s: %d genes, %d sets, %d at padj<0.05 -> %s\n",
            basename(a[1]), length(stats), length(gmt),
            sum(res$padj < 0.05, na.rm = TRUE), basename(a[3])))
