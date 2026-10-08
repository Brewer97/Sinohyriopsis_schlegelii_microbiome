rm(list=ls())

# working directory: the folder that contains the input files
setwd(".")

library(tidyverse)
library(ANCOMBC)
library(phyloseq)

# significance thresholds, following the criteria given in the Methods
P_CUTOFF   <- 0.05   # nominal p value
Q_CUTOFF   <- 0.05   # FDR-adjusted q value
LFC_CUTOFF <- 1      # absolute log2 fold change

# read data
otu <- read.csv("data/otu.csv", row.names = 1)
group <- read.csv("data/group.csv", row.names = 1)
tax <- read.csv("data/tax.csv", row.names = 1)

# build phyloseq object
ps <- phyloseq(
  otu_table(otu, taxa_are_rows = TRUE),
  sample_data(group),
  tax_table(as.matrix(tax))
)

# ANCOMBC2 at genus level
genus_ancombc2_analysis <- function(ps_object, group_var = "Part") {
  cat("Running ANCOMBC2 differential analysis at genus level...\n")
  
  ancombc2_results <- ancombc2(
    data = ps_object,
    assay_name = "counts",
    tax_level = "Genus",
    fix_formula = group_var,
    group = group_var,
    p_adj_method = "fdr",
    prv_cut = 0.10,
    lib_cut = 0,
    struc_zero = TRUE,
    neg_lb = FALSE,
    alpha = 0.05,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE
  )
  
  cat("ANCOMBC2 finished.\n")
  return(ancombc2_results)
}

# run the analysis
ancombc2_results <- genus_ancombc2_analysis(ps, "Part")

# extract genus-level results
extract_genus_results_corrected <- function(ancombc2_obj) {
  cat("=== extracting ANCOMBC2 results ===\n")
  
  if (is.null(ancombc2_obj$res)) {
    stop("ANCOMBC2 did not return the main result table")
  }
  
  res <- ancombc2_obj$res
  cat("result dimensions:", dim(res), "\n")
  
  # result table carrying every column we need
  genus_results <- data.frame(
    Taxon = res$taxon,
    stringsAsFactors = FALSE
  )
  
  # pull the genus name out of the taxon string
  genus_results$Genus <- sapply(strsplit(genus_results$Taxon, "\\|"), function(x) {
    if (length(x) >= 8) {
      return(x[8])  # genus is rank 8
    } else {
      return(tail(x, 1))  # otherwise take the last element
    }
  })
  
  # columns of the PartStom contrast
  genus_results$lfc <- res$lfc_PartStom
  genus_results$p_value <- res$p_PartStom
  genus_results$q_value <- res$q_PartStom
  genus_results$diff <- res$diff_PartStom
  genus_results$passed_ss <- res$passed_ss_PartStom
  genus_results$se <- res$se_PartStom
  genus_results$W <- res$W_PartStom
  
  # significance flags
  genus_results$Significant_p <- genus_results$p_value < P_CUTOFF
  genus_results$Significant_q <- genus_results$q_value < Q_CUTOFF
  
  # sort by p value
  genus_results <- genus_results[order(genus_results$p_value), ]
  
  return(genus_results)
}

# main
tryCatch({
  # genus results
  genus_results <- extract_genus_results_corrected(ancombc2_results)
  
  cat("\n=== summary ===\n")
  cat("total genera:", nrow(genus_results), "\n")
  
  # number of reliable taxa
  reliable_count <- sum(genus_results$passed_ss == TRUE, na.rm = TRUE)
  cat("reliable taxa after sensitivity analysis:", reliable_count, "\n")
  
  # number of significant taxa among the reliable ones
  reliable_sig <- genus_results[genus_results$passed_ss == TRUE & genus_results$q_value < Q_CUTOFF, ]
  cat("reliable and significant taxa (q<0.05):", nrow(reliable_sig), "\n")
  
  # save the reliable taxa
  reliable_results <- genus_results[genus_results$passed_ss == TRUE, ]
  write.csv(reliable_results, "sp_diff/genus_ancombc2_reliable_results.csv", row.names = FALSE)
  cat("genus_ancombc2_reliable_results.csv saved\n")
  
  # taxa retained as significantly and robustly enriched, i.e. the criteria
  # reported in the Methods: FDR < Q_CUTOFF, |log2 fold change| > LFC_CUTOFF,
  # and passed the built-in sensitivity analysis
  enriched_results <- genus_results[genus_results$passed_ss == TRUE &
                                      genus_results$q_value < Q_CUTOFF &
                                      abs(genus_results$lfc) > LFC_CUTOFF, ]
  write.csv(enriched_results, "sp_diff/genus_ancombc2_enriched_results.csv", row.names = FALSE)
  cat("genus_ancombc2_enriched_results.csv saved:", nrow(enriched_results), "genera\n")
  
}, error = function(e) {
  cat("Error:", e$message, "\n")
})
