#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(fgsea); library(data.table); library(jsonlite)})
d <- fread("results/figS1F_gsea_stats_tss.tsv")
gmt <- gmtPathways("results/GO_BP_2025.gmt")
gmt <- lapply(gmt, function(g) unique(sub(",.*", "", trimws(g))))
set.seed(1)
out <- list()
for (kind in c("crush", "eigen")) for (rep in c("1", "2")) {
  col <- sprintf("%s_prom_d%s", kind, rep)
  s <- d[[col]]; names(s) <- d$gene
  s <- s[is.finite(s)]; s <- s[!duplicated(names(s))]
  r <- suppressWarnings(fgseaSimple(gmt, s, nperm = 1000, minSize = 15,
                                    maxSize = 500, scoreType = "std", nproc = 1))
  r <- r[order(pval)]
  ids <- regmatches(r$pathway, regexpr("GO:[0-9]+", r$pathway))
  out[[paste0(kind, rep)]] <- head(ids[!is.na(ids)], 200)
  cat(sprintf("  %s rep%s: %d genes ranked, %d terms scored\n",
              kind, rep, length(s), nrow(r)))
}
write_json(out, "results/top200_goids_prom.json", auto_unbox = FALSE)
cat("wrote results/top200_goids_prom.json\n")
