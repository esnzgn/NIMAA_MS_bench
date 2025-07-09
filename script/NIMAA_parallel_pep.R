rm(list = ls())

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

colnames()

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

  # Get only matching columns
  peptide_cols <- grep("number_of_peptides_[AB]_", colnames(tmptbl), value = TRUE)

  if (length(peptide_cols) == 0) {
    warning(paste("No peptide columns found in", nm))
    return(NULL)
  }

  new_df <- tmptbl %>% select(all_of(peptide_cols))

  result <- extractSubMatrix(
    as.matrix(new_df),
    shape = "Rectangular_element_max",
    col.vars = "samples",
    row.vars = "Proteins"
  )

  assign(nm, result[[1]])
  # list(name = nm, result = result)
}

# Stop cluster after execution
stopCluster(cl)

# Combine results into a named list
subMat <- setNames(
  lapply(subMat, function(x) x$result),
  sapply(subMat, function(x) x$name)
)

for (i in seq_along(subMat))
  tmp <- as.matrix(subMat[[i]])
cls1 <- findCluster(tmp,
                    part = 1,
                    method = "all", # all available clustering methods
                    normalization = TRUE, # normalize the input matrix
                    rm_weak_edges = TRUE, # remove the weak edges in graph
                    rm_method = 'delete', # delete the weak edges instead of lowering their weights to 0.
                    threshold = 'median', # Use median of edges' weights as threshold
                    set_remaining_to_1 = TRUE, # set the weights of remaining edges to 1
)











