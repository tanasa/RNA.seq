#!/usr/bin/env Rscript
# FGSEA only. Runs on single DESeq2-style CSV files (Gene_symbol, log2FoldChange).
# Ranks genes by log2FoldChange (descending). Writes to NEW folders: <BASE>_fgsea/

library(fgsea)
library(readr)
library(dplyr)
library(ggplot2)

# ============ CONFIGURATION ============
# Run from GSEA_ORA folder (script and GMT files in same directory).
# FGSEA input: single CSV per contrast (Ast_Rgl.csv, Gli_IPC.csv, Imm_CGE.csv, ...)
FGSEA_FILES <- c("Ast_Rgl.csv", "Gli_IPC.csv", "Imm_CGE.csv", "Imm_LGE.csv", "Imm_MGE.csv", "Mcg.csv", "OPC.csv")

# GMT files in same folder (GSEA_ORA)
GMT_FILES <- c(
  "pathways_hallmark.gmt",
  "pathways_kegg.gmt",
  "pathways_reactome.gmt",
  "pathways_gobp.gmt",
  "pathways_wikipathway.gmt"
)

# ============ HELPER ============
read_deseq2_csv <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  d <- read_csv(path, col_types = cols(), show_col_types = FALSE)
  if (!"Gene_symbol" %in% colnames(d)) {
    stop("Expected column 'Gene_symbol' in ", path, ". Columns: ", paste(colnames(d), collapse = ", "))
  }
  if (!"log2FoldChange" %in% colnames(d)) {
    stop("Expected column 'log2FoldChange' in ", path)
  }
  d
}

# Load pathway gene sets from GMT files (same folder as script)
pathways <- list()
for (gmt in GMT_FILES) {
  if (file.exists(gmt)) {
    pw <- gmtPathways(gmt)
    # Prefix names to avoid clashes across GMT files
    prefix <- sub("^pathways_|\\.gmt$", "", gmt)
    names(pw) <- paste0(prefix, "::", names(pw))
    pathways <- c(pathways, pw)
    cat("Loaded", length(pw), "pathways from", gmt, "\n")
  } else {
    warning("GMT file not found: ", gmt)
  }
}
cat("Total pathway sets:", length(pathways), "\n")

# ============ LOOP OVER EACH FILE ============
for (csv_file in FGSEA_FILES) {
  if (!file.exists(csv_file)) {
    warning("Skipping ", csv_file, ": file not found")
    next
  }

  # Base name for output folder (e.g. Imm_CGE.csv -> Imm_CGE)
  BASE_NAME <- sub("\\.csv$", "", csv_file)
  # New folder for FGSEA results (do not overwrite existing _enrichR_gsea)
  RESULTS_DIR <- paste0(BASE_NAME, "_fgsea")
  dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

  cat("\n======== ", csv_file, " -> ", RESULTS_DIR, " ========\n", sep = "")

  d <- read_deseq2_csv(csv_file)
  d <- d %>% filter(!is.na(Gene_symbol), Gene_symbol != "", !is.na(log2FoldChange), is.finite(log2FoldChange))
  ranked_vec <- setNames(d$log2FoldChange, d$Gene_symbol)
  ranked_vec <- ranked_vec[!duplicated(names(ranked_vec))]
  ranked_vec <- sort(ranked_vec, decreasing = TRUE)
  cat("Ranked genes:", length(ranked_vec), "\n")

  set.seed(42)
  fgsea_res <- fgsea(pathways = pathways, stats = ranked_vec, minSize = 15, maxSize = 500, eps = 1e-10) %>%
    as.data.frame() %>% arrange(pval)
  write_csv(fgsea_res, file.path(RESULTS_DIR, paste0(BASE_NAME, "_fgsea_all.csv")))

  fgsea_top <- fgsea_res %>% filter(padj < 0.2) %>% head(50)
  if (nrow(fgsea_top) > 0) {
    write_csv(fgsea_top, file.path(RESULTS_DIR, paste0(BASE_NAME, "_fgsea_significant.csv")))
    for (i in seq_len(min(5, nrow(fgsea_top)))) {
      p <- plotEnrichment(pathways[[fgsea_top$pathway[i]]], ranked_vec) +
        labs(title = fgsea_top$pathway[i], subtitle = sprintf("NES=%.2f padj=%.2e", fgsea_top$NES[i], fgsea_top$padj[i]))
      ggsave(file.path(RESULTS_DIR, paste0(BASE_NAME, "_fgsea_", i, "_", gsub("[^A-Za-z0-9_]", "_", fgsea_top$pathway[i]), ".png")), p, width = 8, height = 5, dpi = 150)
    }
    plot_df <- fgsea_top %>% head(20) %>%
      mutate(pathway_short = ifelse(nchar(pathway) > 45, paste0(substr(pathway, 1, 42), "..."), pathway))
    p_bar <- ggplot(plot_df, aes(x = NES, y = reorder(pathway_short, NES), fill = NES > 0)) +
      geom_col() + scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"), guide = "none") +
      labs(title = paste("GSEA -", BASE_NAME), x = "NES", y = "") + theme_minimal() + theme(axis.text.y = element_text(size = 9))
    ggsave(file.path(RESULTS_DIR, paste0(BASE_NAME, "_fgsea_barplot.png")), p_bar, width = 10, height = 7, dpi = 300, bg = "white")
  }

  cat(BASE_NAME, "done. Results in", RESULTS_DIR, "\n")
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("FGSEA complete. Output folders: *_fgsea (new; existing _enrichR_gsea kept)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
