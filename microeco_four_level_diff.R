rm(list=ls())

# working directory: the folder that contains the input files
setwd(".")

library(microeco)
library(tidyverse)

# the two-group contrast tested at every taxonomic level
GROUP_VAR  <- "Part"
CONTRAST   <- "PartStom"   # Part has the levels EW and Stom, so lfc = log2(Stom / EW)
P_ADJUST   <- "fdr"
FDR_CUTOFF <- 0.05         # FDR-adjusted p value threshold for significance

# read data
otu   <- read.csv("data/otu.csv",   row.names = 1)
group <- read.csv("data/group.csv", row.names = 1)
tax   <- read.csv("data/tax.csv",   row.names = 1)

dataset <- microtable$new(otu_table = otu,
                          sample_table = group,
                          tax_table = tax)
dataset$tidy_dataset()

# differential abundance is run once per taxonomic level, inside each subset of
# the level above it (Class within Domain, Order within Class, and so on)
hierarchical_diff_analysis <- function(dataset, target_level, split_level, group_var = GROUP_VAR) {
  cat("differential analysis at", target_level, "level, split by", split_level, "...\n")

  split_categories <- unique(dataset$tax_table[[split_level]])
  split_categories <- split_categories[!is.na(split_categories) & split_categories != ""]

  all_results <- list()
  failed <- character(0)

  for (category in split_categories) {
    cat("  ", split_level, ":", category, "...\n")

    taxa_in_category <- rownames(dataset$tax_table)[dataset$tax_table[[split_level]] == category]

    # at least 2 taxa are needed inside the subset
    if (length(taxa_in_category) <= 1) next

    sub_dataset <- dataset$clone()
    sub_dataset$otu_table <- sub_dataset$otu_table[taxa_in_category, , drop = FALSE]
    sub_dataset$tax_table <- sub_dataset$tax_table[taxa_in_category, , drop = FALSE]
    sub_dataset$tidy_dataset()

    if (nrow(sub_dataset$otu_table) <= 1) next

    ok <- tryCatch({
      diff_obj <- trans_diff$new(
        dataset = sub_dataset,
        taxa_level = target_level,
        p_adjust_method = P_ADJUST,
        method = "ancombc2",
        group = group_var
      )

      if (!is.null(diff_obj$res_diff) && nrow(diff_obj$res_diff) > 0) {
        diff_obj$res_diff[[split_level]] <- category
        all_results[[category]] <- diff_obj$res_diff
      }
      TRUE
    }, error = function(e) {
      cat("    skipped", category, ":", conditionMessage(e), "\n")
      FALSE
    })

    if (!ok) failed <- c(failed, category)
  }

  # a subset can fail on zero variances or degenerate data; keep track of which ones
  if (length(failed) > 0) {
    cat("  subsets that could not be tested:", paste(failed, collapse = ", "), "\n")
  }

  if (length(all_results) == 0) return(NULL)

  final_result <- do.call(rbind, all_results)
  rownames(final_result) <- NULL
  final_result
}

class_result  <- hierarchical_diff_analysis(dataset, "Class",  "Domain")
order_result  <- hierarchical_diff_analysis(dataset, "Order",  "Class")
family_result <- hierarchical_diff_analysis(dataset, "Family", "Order")
genus_result  <- hierarchical_diff_analysis(dataset, "Genus",  "Family")

# keep only the contrast of interest, then drop duplicated terminal names
remove_duplicates_enhanced <- function(diff_result) {
  if (!"Factors" %in% colnames(diff_result)) {
    stop("no 'Factors' column in the result table")
  }
  diff_result <- diff_result[diff_result$Factors == CONTRAST, ]

  diff_result$TaxonName <- sapply(strsplit(diff_result$Taxa, "\\|"), function(x) tail(x, 1))
  diff_result <- diff_result[!duplicated(diff_result$TaxonName), ]
  diff_result$TaxonName <- NULL

  diff_result
}

class_result  <- remove_duplicates_enhanced(class_result)
order_result  <- remove_duplicates_enhanced(order_result)
family_result <- remove_duplicates_enhanced(family_result)
genus_result  <- remove_duplicates_enhanced(genus_result)

# attach the four Significance columns to the full taxonomy of the dataset
organize_diff_results_corrected <- function(class_result, order_result, family_result, genus_result, tax_table) {

  tax_full <- tax_table %>%
    select(Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus) %>%
    distinct()

  class_clean <- class_result %>%
    mutate(Class = str_split(Taxa, "\\|", simplify = TRUE)[, 5]) %>%
    select(Class, Class_Significance = Significance) %>%
    distinct(Class, .keep_all = TRUE)

  order_clean <- order_result %>%
    mutate(
      Class = str_split(Taxa, "\\|", simplify = TRUE)[, 5],
      Order = str_split(Taxa, "\\|", simplify = TRUE)[, 6]
    ) %>%
    select(Class, Order, Order_Significance = Significance) %>%
    distinct(Class, Order, .keep_all = TRUE)

  family_clean <- family_result %>%
    mutate(
      Class = str_split(Taxa, "\\|", simplify = TRUE)[, 5],
      Order = str_split(Taxa, "\\|", simplify = TRUE)[, 6],
      Family = str_split(Taxa, "\\|", simplify = TRUE)[, 7]
    ) %>%
    select(Class, Order, Family, Family_Significance = Significance) %>%
    distinct(Class, Order, Family, .keep_all = TRUE)

  genus_clean <- genus_result %>%
    mutate(
      Class = str_split(Taxa, "\\|", simplify = TRUE)[, 5],
      Order = str_split(Taxa, "\\|", simplify = TRUE)[, 6],
      Family = str_split(Taxa, "\\|", simplify = TRUE)[, 7],
      Genus = str_split(Taxa, "\\|", simplify = TRUE)[, 8]
    ) %>%
    select(Class, Order, Family, Genus, Genus_Significance = Significance) %>%
    distinct(Class, Order, Family, Genus, .keep_all = TRUE)

  # join from the most detailed level upwards
  result <- tax_full %>%
    left_join(genus_clean, by = c("Class", "Order", "Family", "Genus")) %>%
    left_join(family_clean, by = c("Class", "Order", "Family")) %>%
    left_join(order_clean, by = c("Class", "Order")) %>%
    left_join(class_clean, by = "Class")

  result %>%
    select(Domain, Supergroup, Division, Subdivision, Class, Class_Significance,
           Order, Order_Significance,
           Family, Family_Significance,
           Genus, Genus_Significance)
}

final_result_corrected <- organize_diff_results_corrected(class_result, order_result,
                                                          family_result, genus_result, tax)

result <- final_result_corrected[, -c(2:4)]

# a taxon is kept if at least one of its four levels is significant
has_significance <- rowSums(
  !is.na(result[, c("Class_Significance", "Order_Significance",
                    "Family_Significance", "Genus_Significance")]) &
    result[, c("Class_Significance", "Order_Significance",
               "Family_Significance", "Genus_Significance")] != ""
) > 0

result_filtered <- result[has_significance, ]

# count on the FDR-adjusted p value; microeco's Significance marker uses the
# same threshold, so the two give identical counts
count_significant <- function(diff_result) sum(diff_result$P.adj < FDR_CUTOFF, na.rm = TRUE)

cat("=== differential analysis summary ===\n")
cat("significant taxa at class level: ", count_significant(class_result), "\n")
cat("significant taxa at order level: ", count_significant(order_result), "\n")
cat("significant taxa at family level:", count_significant(family_result), "\n")
cat("significant taxa at genus level: ", count_significant(genus_result), "\n")
cat("taxa kept in the filtered table: ", nrow(result_filtered), "\n")

# save
write.csv(class_result,  "sp_diff/class_diff_results.csv")
write.csv(order_result,  "sp_diff/order_diff_results.csv")
write.csv(family_result, "sp_diff/family_diff_results.csv")
write.csv(genus_result,  "sp_diff/genus_diff_results.csv")
write.csv(result_filtered, "sp_diff/final_filtered_results.csv")

cat("differential abundance tables saved\n")
