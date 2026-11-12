#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(fgsea); library(data.table); library(GOSemSim); library(GO.db)})

TOPN <- as.integer(Sys.getenv("TOPN", "100"))
TOPC <- as.integer(Sys.getenv("TOPC", "25"))
d <- fread("results/figS1F_gsea_stats_tss.tsv")
gmt <- gmtPathways("results/GO_BP_2025.gmt")
gmt <- lapply(gmt, function(g) unique(sub(",.*", "", trimws(g))))
set.seed(1)

run <- function(col) {
  s <- d[[col]]; names(s) <- d$gene
  s <- s[is.finite(s)]; s <- s[!duplicated(names(s))]
  r <- suppressWarnings(fgseaSimple(gmt, s, nperm = 1000, minSize = 15,
                                    maxSize = 500, scoreType = "std", nproc = 1))
  r <- r[order(pval)]; r[, rank := .I]; r
}
cr <- run("crush_prom_delta")
invisible(run("crush_prom_d1")); invisible(run("crush_prom_d2"))
eg <- run("eigen_prom_delta")

m <- cr[1:TOPC, .(pathway, CRUSH_NES = NES, CRUSH_p = pval, size,
                  CRUSH_rank = rank)]
m <- merge(m, eg[, .(pathway, eigen_NES = NES, eigen_rank = rank)],
           by = "pathway", all.x = TRUE)
m <- m[order(CRUSH_rank)]

hs <- godata('org.Hs.eg.db', ont = "BP", computeIC = FALSE)
gid <- function(x) regmatches(x, regexpr("GO:[0-9]+", x))
etop <- head(eg, TOPN)
eids <- gid(etop$pathway); eids <- eids[!is.na(eids)]
S <- mgoSim(gid(m$pathway), eids, semData = hs, measure = "Wang", combine = NULL)
S[!is.finite(S)] <- 0
bi <- apply(S, 1, which.max)
m[, wang_best  := round(apply(S, 1, max), 3)]
m[, wang_match := sub("\\s*\\(GO:.*", "", etop$pathway[bi])]
m[, term := sub("\\s*\\(GO:.*", "", pathway)]

out <- "results/figS1H_crush_vs_eigen_averaged.tsv"
fwrite(m[, .(CRUSH_rank, term, CRUSH_NES, CRUSH_p, size,
             eigen_NES, eigen_rank, wang_best, wang_match)],
       out, sep = "\t")

cat(sprintf("%4s %-46s %7s %6s | %7s %6s | %5s %-32s\n",
            "rk", "term", "NES", "n", "eigNES", "eigRk", "wang",
            "closest eigenvector term"))
for (i in seq_len(nrow(m)))
  cat(sprintf("%4d %-46s %7.2f %6d | %7.2f %6d | %5.2f %-32s\n",
              m$CRUSH_rank[i], substr(m$term[i], 1, 46), m$CRUSH_NES[i],
              m$size[i], m$eigen_NES[i], m$eigen_rank[i], m$wang_best[i],
              substr(m$wang_match[i], 1, 32)))
cat(sprintf("\n  median eigenvector rank %.0f of %d terms\n",
            median(m$eigen_rank), nrow(m)))
cat(sprintf("wrote %s\n", out))
