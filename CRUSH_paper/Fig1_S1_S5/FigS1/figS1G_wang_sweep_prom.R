#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(GOSemSim); library(jsonlite); library(GO.db)})

d  <- fromJSON("results/top200_goids_prom.json")
hs <- godata('org.Hs.eg.db', ont = "BP", computeIC = FALSE)
NS <- as.integer(strsplit(Sys.getenv("NS", "10,50,100"), ",")[[1]])
NP <- as.integer(Sys.getenv("NPERM", "60"))

sim <- function(a, b) mgoSim(a, b, semData = hs, measure = "Wang", combine = "BMA")
univ <- keys(GO.db, keytype = "GOID"); univ <- univ[Ontology(GOTERM[univ]) == "BP"]

res <- list()
cat(sprintf("%5s %-26s %8s %8s %8s %8s\n",
            "N", "comparison", "Wang", "null", "nullSD", "z"))
for (n in NS) {
  set.seed(1)
  nl <- replicate(NP, sim(sample(univ, n), sample(univ, n)))
  nl <- nl[is.finite(nl)]
  mu <- mean(nl); sdv <- sd(nl)
  pairs <- list(cross_rep1 = c("crush1","eigen1"), cross_rep2 = c("crush2","eigen2"),
                crush_reps = c("crush1","crush2"), eigen_reps = c("eigen1","eigen2"))
  for (nm in names(pairs)) {
    p <- pairs[[nm]]
    v <- sim(d[[p[1]]][1:n], d[[p[2]]][1:n])
    res[[sprintf("%d|%s", n, nm)]] <- list(N = n, comparison = nm, wang = v,
      null_mean = mu, null_sd = sdv, z = (v - mu) / sdv,
      null_q95 = as.numeric(quantile(nl, 0.95)))
    cat(sprintf("%5d %-26s %8.3f %8.3f %8.3f %8.1f\n", n, nm, v, mu, sdv, (v-mu)/sdv))
  }
  res[[sprintf("%d|null_draws", n)]] <- nl
}
write_json(res, "results/figS1G_wang_sweep_prom.json", auto_unbox = TRUE, digits = 6)
cat("\nwrote results/figS1G_wang_sweep_prom.json\n")
