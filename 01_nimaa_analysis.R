# ==== Load Packages ====
library(tidyverse)
library(NIMAA)
library(igraph)

# ==== Load Data ====
data_dir <- "../WOMBAT-P_Processed/ProteoBenchDDA/0.9.11/"
files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

# Select peptide–protein relationship from multiple tools (example)
file_paths <- list(
  MaxQuant  =  files[grep("stand_pep_quant_mergedmaxquant.csv", files)],
  Proline   =  files[grep("stand_pep_quant_mergedproline.csv", files)],
  CompOmics =  files[grep("stand_pep_quant_mergedcompomics.csv", files)],
  TPP       =  files[grep("stand_pep_quant_mergedtpp.csv", files)]
)

datasets <- map(file_paths, read_csv)

# ==== Convert to Bipartite Graphs ====
nimaa_inputs <- map(datasets, ~ .x %>%
                      select(modified_peptide, protein_group) %>%
                      distinct() %>%
                      drop_na()
)

# ==== Create Bipartite Graphs and Analyze ====
results <- map2(nimaa_inputs, names(nimaa_inputs), function(df, tool) {
  cat("\nAnalyzing", tool, "\n")

  # Plot input matrix
  plotInput(df, plot = TRUE, save = FALSE)

  # Bipartite graph object
  bip <- plotBipartite(df, plot = TRUE, save = FALSE)

  # Submatrix extraction (non-missing)
  submat <- extractSubMatrix(df, method = "row")  # also try method = "column"

  # Cluster projected network
  clust_res <- findCluster(submat, plot = TRUE)

  # Store results
  list(
    bipartite = bip,
    submatrix = submat,
    clustering = clust_res
  )
})

# ==== Optional: Impute Missing Values & Compare ====
impute_results <- map(nimaa_inputs, function(df) {
  imputeMissingValue(df, method = c("mean", "CA", "als", "svd"), plot = TRUE)
})

# ==== Save visualizations ====
dir.create("results", showWarnings = FALSE)
# (You can export plots using ggsave or pdf() if needed)

# ==== Example Summary of Cluster Quality ====
map(results, ~ .x$clustering$metrics) %>%
  bind_rows(.id = "Tool") %>%
  arrange(desc(Modularity)) %>%
  print()
