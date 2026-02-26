# ============================================================================
# SCRIPT: DOWNLOAD PATHWAY FILES FOR GSEA
# ============================================================================
# Downloads pathway gene sets from MSigDB
# ============================================================================

library(msigdbr)
library(fgsea)

cat("Downloading pathway gene sets from MSigDB...\n\n")

# Function to fetch pathways using msigdbr
fetch_pathways <- function(collection = NULL, subcollection = NULL) {
  if (is.null(collection)) {
    stop("Collection must be specified")
  }
  
  # Fetch using msigdbr (using new API with collection/subcollection)
  if (collection == "H") {
    # Hallmark pathways
    pathways <- msigdbr(species = "Homo sapiens", collection = "H")
    pathway_list <- split(pathways$gene_symbol, pathways$gs_name)
    return(pathway_list)
  } else if (collection == "C2") {
    # C2: Curated gene sets
    if (!is.null(subcollection)) {
      if (subcollection == "CP:REACTOME") {
        pathways <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
      } else if (subcollection == "CP:KEGG_MEDICUS") {
        # Try KEGG_MEDICUS first, then fall back to KEGG
        tryCatch({
          pathways <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS")
        }, error = function(e) {
          cat("   Note: KEGG_MEDICUS not found, trying CP:KEGG\n")
          pathways <<- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG")
        })
      } else if (subcollection == "CP:WIKIPATHWAYS") {
        pathways <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS")
      } else {
        pathways <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = subcollection)
      }
    } else {
      pathways <- msigdbr(species = "Homo sapiens", collection = "C2")
    }
    pathway_list <- split(pathways$gene_symbol, pathways$gs_name)
    return(pathway_list)
  } else if (collection == "C5") {
    # C5: GO gene sets
    if (!is.null(subcollection)) {
      if (subcollection == "GO:BP") {
        pathways <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
      } else {
        pathways <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = subcollection)
      }
    } else {
      pathways <- msigdbr(species = "Homo sapiens", collection = "C5")
    }
    pathway_list <- split(pathways$gene_symbol, pathways$gs_name)
    return(pathway_list)
  } else {
    pathways <- msigdbr(species = "Homo sapiens", collection = collection)
    pathway_list <- split(pathways$gene_symbol, pathways$gs_name)
    return(pathway_list)
  }
}

# Create output directory
dir.create("pathways", showWarnings = FALSE, recursive = TRUE)

# Download pathways
cat("1. Downloading Hallmark pathways (H)...\n")
pathways_hallmark <- fetch_pathways("H")
cat("   Downloaded", length(pathways_hallmark), "Hallmark pathways\n")
saveRDS(pathways_hallmark, "pathways/pathways_hallmark.rds")
cat("   Saved to: pathways/pathways_hallmark.rds\n\n")

cat("2. Downloading Reactome pathways (C2, CP:REACTOME)...\n")
pathways_reactome <- fetch_pathways("C2", "CP:REACTOME")
cat("   Downloaded", length(pathways_reactome), "Reactome pathways\n")
saveRDS(pathways_reactome, "pathways/pathways_reactome.rds")
cat("   Saved to: pathways/pathways_reactome.rds\n\n")

cat("3. Downloading KEGG pathways (C2, CP:KEGG_MEDICUS)...\n")
# Note: KEGG_MEDICUS might not be available, using CP:KEGG instead
pathways_kegg <- tryCatch({
  fetch_pathways("C2", "CP:KEGG_MEDICUS")
}, error = function(e) {
  cat("   Note: KEGG_MEDICUS not found, using CP:KEGG instead\n")
  fetch_pathways("C2", "CP:KEGG")
})
cat("   Downloaded", length(pathways_kegg), "KEGG pathways\n")
saveRDS(pathways_kegg, "pathways/pathways_kegg.rds")
cat("   Saved to: pathways/pathways_kegg.rds\n\n")

cat("4. Downloading GO Biological Process pathways (C5, GO:BP)...\n")
pathways_gobp <- fetch_pathways("C5", "GO:BP")
cat("   Downloaded", length(pathways_gobp), "GO:BP pathways\n")
saveRDS(pathways_gobp, "pathways/pathways_gobp.rds")
cat("   Saved to: pathways/pathways_gobp.rds\n\n")

cat("5. Downloading WikiPathways (C2, CP:WIKIPATHWAYS)...\n")
pathways_wikipathway <- fetch_pathways("C2", "CP:WIKIPATHWAYS")
cat("   Downloaded", length(pathways_wikipathway), "WikiPathways\n")
saveRDS(pathways_wikipathway, "pathways/pathways_wikipathway.rds")
cat("   Saved to: pathways/pathways_wikipathway.rds\n\n")

# Also save as .gmt format for compatibility
cat("Converting to .gmt format...\n")
write_gmt <- function(pathway_list, filename) {
  con <- file(filename, "w")
  for (pathway_name in names(pathway_list)) {
    genes <- pathway_list[[pathway_name]]
    line <- paste(pathway_name, "NA", paste(genes, collapse = "\t"), sep = "\t")
    writeLines(line, con)
  }
  close(con)
}

write_gmt(pathways_hallmark, "pathways/pathways_hallmark.gmt")
write_gmt(pathways_reactome, "pathways/pathways_reactome.gmt")
write_gmt(pathways_kegg, "pathways/pathways_kegg.gmt")
write_gmt(pathways_gobp, "pathways/pathways_gobp.gmt")
write_gmt(pathways_wikipathway, "pathways/pathways_wikipathway.gmt")

cat("   Saved .gmt files to pathways/ directory\n\n")

cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("DOWNLOAD COMPLETE!\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")
cat("Pathways saved as:\n")
cat("  - RDS files (for R): pathways/*.rds\n")
cat("  - GMT files (for GSEA): pathways/*.gmt\n\n")
cat("Pathway counts:\n")
cat("  - Hallmark:", length(pathways_hallmark), "\n")
cat("  - Reactome:", length(pathways_reactome), "\n")
cat("  - KEGG:", length(pathways_kegg), "\n")
cat("  - GO:BP:", length(pathways_gobp), "\n")
cat("  - WikiPathways:", length(pathways_wikipathway), "\n")

