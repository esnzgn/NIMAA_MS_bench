# # Make sure devtools is installed
# install.packages("devtools")
#
# # Install NIMAA from GitHub (latest features)
# devtools::install_github("jafarilab/NIMAA")

# ==== Load Packages ====
library(tidyverse)
library(dplyr)
library(NIMAA)
library(igraph)
library(purrr)
library(readr)

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

# new POV ====
# n*n of peptide to samples
for (i in 1:length(datasets)){
  cat(" ------> ", names(datasets[i]))
  cat("\n colnames ", colnames(datasets[[i]]))
  print(head(datasets[[i]]))
}

# ==== Convert to Bipartite Graphs ====
nimaa_inputs <- map(datasets, ~ .x %>%
                      select(modified_peptide, protein_group, abundance_A_1, abundance_A_2,
                             abundance_A_3, abundance_B_1, abundance_B_2, abundance_B_3) %>%
                      distinct() %>%
                      drop_na()
)

for (i in 1:length(nimaa_inputs)){
  cat(" ------> ", names(nimaa_inputs[i]))
  cat("\n colnames ", colnames(nimaa_inputs[[i]]))
  print(head(nimaa_inputs[[i]]))
}

# --- Define normalization function (Z-score per row or globally) ---
z_score <- function(x) {
  if (all(is.na(x))) return(rep(NA, length(x)))
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# --- Apply to all tools ---
nimaa_long_inputs <- map(nimaa_inputs, function(df) {
  # Select numeric abundance columns
  abund_cols <- grep("^abundance_", colnames(df), value = TRUE)

  # Normalize across all samples (global Z-score per row)
  df_norm <- df %>%
    mutate(across(all_of(abund_cols), z_score)) %>%
    select(modified_peptide, all_of(abund_cols)) %>%
    pivot_longer(
      cols = all_of(abund_cols),
      names_to = "sample",
      values_to = "z_abundance"
    ) %>%
    mutate(sample = gsub("abundance_", "", sample))

  return(df_norm)
})

for (i in 1:length(nimaa_long_inputs)){
  print(head())
}

inc_matrix <- NIMAA::plotIncMatrix(x = data.frame(nimaa_long_inputs[[1]]),
                                   index_nominal = c(1, 2),
                                   index_numeric = "z_abundance"
)


inc_matrix <- createIncidenceMatrix(
  df_long,
  part1 = "modified_peptide",
  part2 = "sample",
  value = "z_abundance"
)

plotIncMatrix(inc_matrix)




# Compare missingness structure
missing_summary <- map_df(nimaa_analysis, function(res) {
  m <- res$inc$matrix
  tibble(
    Tool = res$tool,
    MissingPercent = mean(is.na(m)) * 100,
    Rows = nrow(m),
    Cols = ncol(m)
  )
})

print(missing_summary)

