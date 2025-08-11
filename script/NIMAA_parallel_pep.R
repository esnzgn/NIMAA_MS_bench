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
library(igraph)

#----
data_dir <- "../WOMBAT-P_Processed/ProteoBenchDDA/0.9.11/"

# List all files
all_files <- list.files(data_dir, full.names = TRUE)

# Split by type
csv_files <- all_files[grepl("\\.(csv)$", all_files)]

# Read files
csv_data <- setNames(lapply(csv_files, read_csv), basename(csv_files))

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
# # clusterEvalQ(cl, {
# #   user_lib <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library", paste0(R.version$major, ".", R.version$minor))
# #   .libPaths(c(user_lib, .libPaths()))
# # })
#
# # Register the parallel backend
# registerDoParallel(cl)

# subMat <- foreach(i = seqq, .packages = c("tidyverse", "NIMAA")) %dopar% {
for(i in seqq){
  nm <- names(csv_data[i])
  cat("---------------->>  prot csv file: ", nm)
  tmptbl <- csv_data[[i]]
  peptide_cols <-  ("number_of_peptides_[AB]_", colnames(tmptbl), value = TRUE)

  # Get only matching columns

  if (length(peptide_cols) == 0) {
    warning(paste("No peptide columns found in", nm))
    return(NULL)
  }

  new_df <- tmptbl %>% select(all_of(peptide_cols))

  if (grepl("stand_prot_quant_mergedmaxquant.csv",nm)){
      for (col_i in 1:ncol(new_df)){
        zero_ids <- which(new_df[, col_i] == 0)
        new_df[zero_ids, col_i] <- NA
      }
  }

  result <- extractSubMatrix(
    as.matrix(new_df),
    shape = "Rectangular_element_max",
    col.vars = "samples",
    row.vars = "Proteins"
  )

  assign(nm, result[[1]])
  # list(name = nm, result = result)
}

Rect_max_prot <- tibble(dataset = c("compomics", "maxquant", "proline", "tpp"),
                   info_lvl = "Protein",
                   shape = "Rectangular_element_max",
                   nestedness_rect_max = c(23.27587, 3.947907, 16.72077, 16.20169),
                   nrow_rect_max = c(2612, 4500, 3547, 3976),
                   ncol_rect_max = 6
                   )

non_missing_perc_vs_ttl <- NULL
for (i in 1:4){
  nrrow_ttl_iden <- nrow(csv_data[[i+8]])
  non_missing_perc_vs_ttl <- c(non_missing_perc_vs_ttl,
                              (Rect_max_prot$nrow_rect_max[i] *
                               Rect_max_prot$ncol_rect_max[i]) / (nrrow_ttl_iden * 6))
}

Rect_max_prot <- cbind(Rect_max_prot, non_missing_perc_vs_ttl)

relative_compl <- NULL
intial_nrows <- NULL
for(i in seqq){
  nm <- names(csv_data[i])
  cat("---------------->>  prot csv file: ", nm)
  tmptbl <- csv_data[[i]]
  peptide_cols <- grep("number_of_peptides_[AB]_", colnames(tmptbl), value = TRUE)

  # Get only matching columns

  if (length(peptide_cols) == 0) {
    warning(paste("No peptide columns found in", nm))
    return(NULL)
  }

  new_df <- tmptbl %>% select(all_of(peptide_cols))

  if (grepl("stand_prot_quant_mergedmaxquant.csv",nm)){
    for (col_i in 1:ncol(new_df)){
      zero_ids <- which(new_df[, col_i] == 0)
      new_df[zero_ids, col_i] <- NA
    }
  }

  tmp_nm <- sum(is.na(new_df)) / (nrow(new_df)[1] * ncol(new_df)[1])

  relative_compl <- c(relative_compl, (1 - tmp_nm))
  intial_nrows <- c(intial_nrows, nrow(new_df))
}

Rect_max_prot$rel_non_missing <- relative_compl
Rect_max_prot$intial_nrow <- intial_nrows
seqq_sub_mat <- ls(pattern = "stand_prot_quant_")

# for(i in seqq){
#   nm <- names(csv_data[i])
#   tmptbl <- csv_data[[i]] %>%
#     data_frame()
#   peptide_cols <- grep("number_of_peptides_[AB]_", colnames(tmptbl), value = TRUE)
#
#
#   # Get only matching columns
#
#   if (length(peptide_cols) == 0) {
#     warning(paste("No peptide columns found in", nm))
#     return(NULL)
#   }
#
#   cat("---------------->>  prot csv file: ", nm, " --- dim:  ",dim((tmptbl[, peptide_cols])), " \n")
#
#   bin_tmptbl <- ifelse(is.na((tmptbl[, peptide_cols])), yes = 0, no = 1)
#
#   # calling a func to cal the clustering
#   tmp <- hclust(bin_tmptbl, method = "ward.D2")
#
#   # sum of the no of clusters to be added to the last col of the rect_max featur
#
# }

# 4 mansucript
## adding the vis of rectmax rowvis and table of whiteboard
library(ggpubr)
library(dplyr)

Rect_max_prot_clean <- Rect_max_prot %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  rename_with(~ gsub("_", " ", .x))  # Replace underscores with spaces for cleaner labels

Rect_max_prot_clean
p_table <- ggtexttable(Rect_max_prot_clean,
                       rows = NULL,  # remove row numbers
                       theme = ttheme("light"))  # or try "mOrange", "classic", "minimal"

# Show table
print(p_table)
p_table
# Save as image (e.g., PNG for publication)
ggsave("Rectangular_Protein_Table.png", p_table, width = 1500, height = 500, units = "px", dpi = 300)

# Or as PDF
ggsave("Rectangular_Protein_Table.pdf", p_table, width = 10, height = 4)

## adding the text for those above







