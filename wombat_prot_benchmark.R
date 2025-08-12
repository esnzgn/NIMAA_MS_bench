# rm(list = ls())

# Define a path in your home directory (e.g., ~/R/4.4)
# user_lib <- file.path(Sys.getenv("HOME"), "R", paste0("x86_64-pc-linux-gnu-library"), paste0(R.version$major, ".", R.version$minor))

# Create the directory if it doesn't exist
# dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)

# Add it to my library paths
# .libPaths(user_lib[[1]])

# Confirm
# .libPaths()

# install.packages("NIMAA")  # or BiocManager::install("NIMAA") if it's from Bioconductor

######
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(jsonlite)
library(readr)
library(purrr)
library(tibble)
library(tidyverse)
library(ggplot2)
library(NIMAA)
library(janitor)
library(doParallel)
library(foreach)

#----
data_dir <- "../WOMBAT-P_Processed/ProteoBenchDDA/0.9.11/"

# List all files
all_files <- list.files(data_dir, full.names = TRUE)

# Split by type
csv_files <- all_files[grepl("\\.(csv)$", all_files)]
json_files <- all_files[grepl("\\.(json)$", all_files)]

# Read files
csv_data <- setNames(lapply(csv_files, read_csv), basename(csv_files))
json_data <- setNames(lapply(json_files, fromJSON), basename(json_files))

# summary stat
# csv_data_log <- csv_data$stand_pep_quant_mergedcompomics.csv %>%
#   mutate(across(starts_with("abundance_"), ~log2(.x + 1)))  # add 1 to avoid log(0)
#
# csv_data_log %>%
#   select(starts_with("abundance_")) %>%
#   summary()
#
# tibble(csv_data$stand_pep_quant_mergedcompomics.csv) %>%
#   select(starts_with("abundance_")) %>%
#   summary()

# convertor func

# pivot long for normalizaiton
for (i in seq_along(csv_data)){
  cat("------> ", i, "\n")
  tmptbl = NULL
  df_long = NULL
  name_file <- sub("\\.csv$", "", sub("^stand_", "", names(csv_data[i])))

  # head(tibble(csv_data[[i]]))
  tmptbl <- csv_data[[i]]

  df_long <- tmptbl %>%
    pivot_longer(
      cols = starts_with("abundance_"),
      names_to = "sample",
      values_to = "abundance"
    )

  p1 <- ggplot(data = df_long, aes(x = sample, y = abundance, fill = sample)) +
    geom_boxplot(alpha = .7) +
    # coord_flip() +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")

  #
  # df_long <- df_long %>%
  #   mutate(across("abundance", ~log2(.x + 1)))

  p2 <- ggplot(data = df_long, aes(x = sample, y = log2(abundance + 1), fill = sample)) +
    geom_boxplot(alpha = .7) +
    # coord_flip() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # legend <- cowplot::get_legend(p1)
  # p1 <- p1 + theme(legened.position = "none")
  p <- cowplot::plot_grid(p1, p2, ncol = 1, labels = gsub("stand_pep_quant_","",names(csv_data[i])))
  print(p)
  ggsave(file.path("./output", paste0(name_file, ".jpg")), plot = p)
}

# func4NiM ----
df_4_NIMAA <- function(not_ready_df, column_to_row, pattern_to_columns){
  # (tmptbl, column_to_row =  "modified_peptide",
  #  pattern_to_columns = c("abundance_A_", "abundance_B_"))
  row_id <- grepl(column_to_row, colnames(not_ready_df))

  ids_rm <- not_ready_df[, which(row_id)][[1]]
  # ids_rm <- c(ids_rm, grep("kit", (not_ready_df[, which(row_id)][[1]]), ignore.case = T))
  # not_ready_df <- not_ready_df[-ids_rm, ]
  # strsplit(not_ready_df[, which(row_id)], split = ",")

  # head(csv_data["stand_pep_quant_mergedcompomics.csv"])
  new_df = NULL

  ids <- grep("abundance_", colnames(not_ready_df))
  # if (length(ids) > 0){
  # not_ready_df <- cbind(ids_rm, not_ready_df[,c(ids)])
  # }
  not_ready_df <- not_ready_df[,c(ids)]
  # if (length(pattern_to_columns) > 1){
  #   cols_ids = NULL
  #   for (ij in 1:length(pattern_to_columns)){
  #     cols_ids <- c(cols_ids, grep(pattern_to_columns[ij], colnames(not_ready_df)))
  #   }
  #   new_df <- not_ready_df[,cols_ids]
  # } else {
  #   new_df <- not_ready_df[,cols_ids]
  # }

  row_names <- ids_rm
  id_nna <- which(!is.na(row_names))
  id_Ndup <- which(!duplicated(row_names))

  id_good <- intersect(id_nna, id_Ndup)

  row_names <- not_ready_df$ids_rm

  new_df <- as.matrix(not_ready_df[id_good, ])
  rownames(new_df) <- row_names
  return(new_df)
}

# analysis ----
# Get the indices of interest
seqq <- grep("stand_prot", names(csv_data))
subMat = NULL

# Set up the cluster
# n_cores <- 8
# cl <- makeCluster(n_cores)
#
# # Manually define the path INSIDE the cluster
# clusterEvalQ(cl, {
#   user_lib <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library", paste0(R.version$major, ".", R.version$minor))
#   .libPaths(c(user_lib, .libPaths()))
# })
#
# # Register the parallel backend
# registerDoParallel(cl)

# subMat <- foreach(i = seqq, .packages = c("tidyverse", "NIMAA")) %dopar% {
idss <- data_frame(NULL)

# Palette for organisms (pick 4 distinct colors)
org_pal <- c(
  "Homo s" = "#1f77b4",
  "Saccha" = "#2ca02c",
  "Escher" = "#d62728",
  "Unknown" = "#7f7f7f"
)

source("./script/extractSubMatrix_complete.R")

for (i in seqq[2:4]){
  nm <- names(csv_data[i])
  tmptbl <- csv_data[[i]]

  new_df <- tmptbl %>%
    select(contains(c("number_of_peptides_A", "number_of_peptides_B")))

  if (grepl("stand_prot_quant_mergedmaxquant.csv",nm)){
    for (col_i in 1:ncol(new_df)){
      zero_ids <- which(new_df[, col_i] == 0)
      new_df[zero_ids, col_i] <- NA
    }
  }

  result <- extractSubMatrix_complete(
    as.matrix(new_df),
    shape = "Rectangular_element_max",
    col.vars = "samples",
    row.vars = "peptides"
  )

  ids_dropped <- result$Rectangular_element_max_dropped %>%
    as.data.frame() %>%
    rownames() %>%
    as.integer()

  # find("extractSubMatrix")
  # body(extractSubMatrix)
  # utils::page("extractSubMatrix")

  ids <-  result$Rectangular_element_max %>%
          as.data.frame() %>%
          rownames() %>%
          as.integer()

  prot_ls <- tmptbl$protein_group[ids]
  prot_ls_dropped <- tmptbl$protein_group[ids_dropped]

  idx_org <- match(prot_ls, result_organism$accession)
  organisms <- substr(result_organism$organism[idx_org], 1, 6)

  idx_org_drp <- match(prot_ls_dropped, result_organism$accession)
  organisms_drp <- substr(result_organism$organism[idx_org_drp], 1, 6)

  forheat <- result$Rectangular_element_max %>%
    as.data.frame() %>%
    cbind(organisms)

  org_vec <- organisms_drp
  org_vec[is.na(org_vec)] <- "Unknown"
  df <- data.frame(row_id = seq_along(org_vec), organism = org_vec)

  # Make sure NA is labeled
  forheat <- forheat %>%
    mutate(organisms = ifelse(is.na(organisms), "Unknown", organisms))

  # Color palette
  org_pal <- c(
    "Homo s" = "#1f77b4",
    "Saccha" = "#2ca02c",
    "Escher" = "#d62728",
    "Unknown" = "#7f7f7f"
  )

  # Add row ID
  df_anno <- forheat %>%
    mutate(row_id = row_number())

  # Plot non_missing part
  p <- ggplot(df_anno, aes(x = 1, y = row_id, fill = organisms)) +
    geom_tile() +
    scale_fill_manual(values = org_pal, name = "Organism") +
    theme_void() +
    theme(
      legend.position = "right",
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
# dropped part
  p_droped <- ggplot(df, aes(x = 1, y = row_id, fill = organism)) +
    geom_tile() +
    scale_fill_manual(
      values = org_pal,
      breaks = c("Homo s","Escher","Saccha","Unknown"),
      labels = c("Human","E. coli","Saccharomyces","Unknown"),
      name = "Organism"
    ) +
    theme_void() +
    theme(legend.position = "right")

  # 5) Save with dynamic height (so tall vectors are readable)
  # h <- max(4, length(org_vec) / 200)  # tweak divisor to control density
  # ggsave(paste0("./output/dropped_organisms_annotation", nm,".png"), p_droped, width = 12, height = h, dpi = 300)
  #
  # # Save image
  # ggsave(paste0("./output/","organism_annotation_", nm,".png"), p, width = 2, height = 12, dpi = 300)
  #
  # write.csv(result$Rectangular_element_max, paste0("./output/", nm, ".csv"))
  cat("##### non missing sub matrix ####")
  print(table(df_anno$organisms))

  cat("##### missing remaining part ####")
  print(table(df$organism))

  cat(" --------> ", nm)
}




# --- prep ---







