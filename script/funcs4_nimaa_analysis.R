df_4_NIMAA <- function(not_ready_df, column_to_row, if_protein = TRUE, pattern_to_columns, No_rows = 1000){
  # df_4_NIMAA(tmptbl, column_to_row = "protein_group", if_protein = T, pattern_to_columns = c("peptides_A", "peptides_B"))
  # df_4_NIMAA(tmptbl, column_to_row = "peptide", pattern_to_columns = c("abundance_A", "abundance_B"))
  # pattern_to_columns <- c("abundance_A", "abundance_B")
  # not_ready_df <- csv_data[[i]]
  # print(names(csv_data[i]))
  # not_ready_df <- csv_data[["stand_pep_quant_mergedcompomics.csv"]]
  # row_id <- grepl("peptide", colnames(csv_data[["stand_pep_quant_mergedcompomics.csv"]]))
  # cols_ids <- grepl("abundance", colnames(csv_data[["stand_pep_quant_mergedcompomics.csv"]]))
  not_ready_df <- not_ready_df[1:No_rows,]
  row_id <- grepl(column_to_row, colnames(not_ready_df))

  if (isTRUE(if_protein)){
    ids_rm <- grep("__", (not_ready_df[, which(row_id)][[1]]))
    ids_rm <- c(ids_rm, grep("kit", (not_ready_df[, which(row_id)][[1]]), ignore.case = T))
    not_ready_df <- not_ready_df[-ids_rm, ]
    # strsplit(not_ready_df[, which(row_id)], split = ",")
  }

  # head(csv_data["stand_pep_quant_mergedcompomics.csv"])
  new_df = NULL

  ids_rm <- grep(("differential_"), colnames(not_ready_df))
  if (length(ids_rm) > 0){
    not_ready_df <- not_ready_df[,-c(ids_rm)]
  }

  if (length(pattern_to_columns) > 1){
    cols_ids = NULL
    for (ij in 1:length(pattern_to_columns)){
      cols_ids <- c(cols_ids, grep(pattern_to_columns[ij], colnames(not_ready_df)))
    }
    new_df <- not_ready_df[,cols_ids]
  } else {
    new_df <- not_ready_df[,cols_ids]
  }

  row_names <- not_ready_df[,which(row_id)]
  id_nna <- which(!is.na(row_names))
  id_Ndup <- which(!duplicated(row_names))

  id_good <- intersect(id_nna, id_Ndup)

  row_names <- row_names[[1]][id_good]

  new_df <- as.data.frame(new_df[id_good, ])
  rownames(new_df) <- row_names
  return(new_df)
}
