# clean the env
rm(list = ls())

# Load required libraries
library(jsonlite)
library(readr)
library(purrr)
library(tibble)

# Set your working directory to the cloned folder
data_dir <- "../WOMBAT-P_Processed/ProteoBenchDDA/0.9.11/"

# List all files in the directory
all_files <- list.files(data_dir, full.names = TRUE)

# Split files by type
csv_files <- all_files[grepl("\\.csv$", all_files)]
json_files <- all_files[grepl("\\.json$", all_files)]

# Read CSV files into a named list
csv_data <- setNames(lapply(csv_files, read_csv), basename(csv_files))

# Read JSON files into a named list
json_data <- setNames(lapply(json_files, fromJSON), basename(json_files))

# Function to print basic metadata
print_metadata <- function(data_list, type = "CSV") {
  cat("\n===== Metadata for", type, "files =====\n")
  for (name in names(data_list)) {
    cat("\n---", name, "---\n")
    df <- data_list[[name]]
    if (is.data.frame(df)) {
      cat("Dimensions:", paste(dim(df), collapse = " x "), "\n")
      cat("Column names:", paste(colnames(df), collapse = ", "), "\n")
      cat("First few row names / IDs:", paste(rownames(df)[1:min(5, nrow(df))], collapse = ", "), "\n")
    } else {
      cat("Not a data frame. Type:", class(df), "\n")
      cat("length is: ", length(df))
      cat("including: ", names(df), "\n")

    }
  }
}

# Print metadata
print_metadata(csv_data, "CSV")
print_metadata(json_data, "JSON")
