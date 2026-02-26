#!/usr/bin/env Rscript
# EnrichR (ORA) only. Reads pos/neg CSV files from current directory.
# Upregulated: *_posLog2FC0.5.csv. Downregulated: *_negLog2FC0.5.csv.
# Output: <BASE_NAME>_enrichR_gsea/ (existing folders kept).

library(enrichR)
library(readr)
library(dplyr)
library(ggplot2)

# ============ CONFIGURATION ============
# CSV files read from current working directory.
# Pos/neg file names: <BASE_NAME>_posLog2FC0.5.csv and _negLog2FC0.5.csv
BASE_NAMES <- c("Imm_CGE", "Imm_LGE", "Imm_MGE", "Mcg", "OPC")

# Expected columns: Gene_symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
read_deseq2_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  d <- read_csv(path, col_types = cols(), show_col_types = FALSE)
  if (!"Gene_symbol" %in% colnames(d) || !"log2FoldChange" %in% colnames(d)) return(NULL)
  list(data = d, gene_col = "Gene_symbol")
}

enrichr_databases <- c(
  "GO_Biological_Process_2025", "KEGG_2021_Human", "Reactome_Pathways_2024",
  "WikiPathways_2024_Human", "MSigDB_Hallmark_2020", "Reactome_2022", "WikiPathways_2023_Human"
)
adj_p_cutoff <- 0.2
top_n <- 50

# ============ LOOP OVER EACH CONTRAST ============
for (BASE_NAME in BASE_NAMES) {
  FILE_POS <- paste0(BASE_NAME, "_posLog2FC0.5.csv")
  FILE_NEG <- paste0(BASE_NAME, "_negLog2FC0.5.csv")
  if (!file.exists(FILE_POS) || !file.exists(FILE_NEG)) {
    warning("Skipping ", BASE_NAME, ": missing ", FILE_POS, " or ", FILE_NEG)
    next
  }

  RESULTS_DIR <- paste0(BASE_NAME, "_enrichR_gsea")
  dir.create(file.path(RESULTS_DIR, "enrichR", "tables"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(RESULTS_DIR, "enrichR", "plots"), showWarnings = FALSE, recursive = TRUE)
  cat("\n======== ", BASE_NAME, " -> ", RESULTS_DIR, " ========\n", sep = "")

  pos <- read_deseq2_csv(FILE_POS)
  pos_df <- pos$data %>% mutate(.gene = .data[[pos$gene_col]]) %>%
    filter(!is.na(.gene), .gene != "", !is.na(log2FoldChange))
  neg <- read_deseq2_csv(FILE_NEG)
  neg_df <- neg$data %>% mutate(.gene = .data[[neg$gene_col]]) %>%
    filter(!is.na(.gene), .gene != "", !is.na(log2FoldChange))

  # Enrichr: upregulated
  genes_up <- unique(pos_df$.gene); genes_up <- genes_up[!is.na(genes_up) & genes_up != ""]
  if (length(genes_up) >= 5) {
    enrichr_up <- enrichr(genes_up, enrichr_databases)
    base_up <- paste0(BASE_NAME, "_upregulated")
    dir.create(file.path(RESULTS_DIR, "enrichR", "tables", base_up), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(RESULTS_DIR, "enrichR", "plots", base_up), showWarnings = FALSE, recursive = TRUE)
    write_csv(data.frame(Gene = genes_up), file.path(RESULTS_DIR, "enrichR", "tables", base_up, paste0(base_up, "_genes.csv")))
    for (database in names(enrichr_up)) {
      res <- enrichr_up[[database]]
      if (nrow(res) > 0) {
        sig <- res %>% filter(Adjusted.P.value < adj_p_cutoff) %>% arrange(Adjusted.P.value) %>% head(top_n)
        if (nrow(sig) > 0) {
          write_csv(sig, file.path(RESULTS_DIR, "enrichR", "tables", base_up, paste0(base_up, "_", database, ".csv")))
          plot_data <- sig %>% head(15) %>% mutate(Term_short = ifelse(nchar(Term) > 50, paste0(substr(Term, 1, 47), "..."), Term), log_p = -log10(Adjusted.P.value))
          p <- ggplot(plot_data, aes(x = log_p, y = reorder(Term_short, log_p))) + geom_col(fill = "firebrick", alpha = 0.7) +
            labs(title = paste("Enrichr (up) -", database), x = "-log10(Adj.P)", y = "") + theme_minimal() + theme(axis.text.y = element_text(size = 9))
          ggsave(file.path(RESULTS_DIR, "enrichR", "plots", base_up, paste0(base_up, "_", database, ".png")), p, width = 10, height = 6, dpi = 300, bg = "white")
        }
      }
    }
  }

  # Enrichr: downregulated
  genes_down <- unique(neg_df$.gene); genes_down <- genes_down[!is.na(genes_down) & genes_down != ""]
  if (length(genes_down) >= 5) {
    enrichr_down <- enrichr(genes_down, enrichr_databases)
    base_down <- paste0(BASE_NAME, "_downregulated")
    dir.create(file.path(RESULTS_DIR, "enrichR", "tables", base_down), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(RESULTS_DIR, "enrichR", "plots", base_down), showWarnings = FALSE, recursive = TRUE)
    write_csv(data.frame(Gene = genes_down), file.path(RESULTS_DIR, "enrichR", "tables", base_down, paste0(base_down, "_genes.csv")))
    for (database in names(enrichr_down)) {
      res <- enrichr_down[[database]]
      if (nrow(res) > 0) {
        sig <- res %>% filter(Adjusted.P.value < adj_p_cutoff) %>% arrange(Adjusted.P.value) %>% head(top_n)
        if (nrow(sig) > 0) {
          write_csv(sig, file.path(RESULTS_DIR, "enrichR", "tables", base_down, paste0(base_down, "_", database, ".csv")))
          plot_data <- sig %>% head(15) %>% mutate(Term_short = ifelse(nchar(Term) > 50, paste0(substr(Term, 1, 47), "..."), Term), log_p = -log10(Adjusted.P.value))
          p <- ggplot(plot_data, aes(x = log_p, y = reorder(Term_short, log_p))) + geom_col(fill = "steelblue", alpha = 0.7) +
            labs(title = paste("Enrichr (down) -", database), x = "-log10(Adj.P)", y = "") + theme_minimal() + theme(axis.text.y = element_text(size = 9))
          ggsave(file.path(RESULTS_DIR, "enrichR", "plots", base_down, paste0(base_down, "_", database, ".png")), p, width = 10, height = 6, dpi = 300, bg = "white")
        }
      }
    }
  }

  cat(BASE_NAME, "done. Results in", RESULTS_DIR, "\n")
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("EnrichR analysis complete. Output folders: *_enrichR_gsea\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
