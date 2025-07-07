df_4_NIMAA <- function(not_ready_df, column_to_row, if_protein = TRUE, pattern_to_columns){
  # df_4_NIMAA(tmptbl, column_to_row = "protein_group", if_protein = T, pattern_to_columns = c("peptides_A", "peptides_B"))
  # df_4_NIMAA(tmptbl, column_to_row = "peptide", pattern_to_columns = c("abundance_A", "abundance_B"))
  # pattern_to_columns <- c("abundance_A", "abundance_B")
  # not_ready_df <- csv_data[[i]]
  # print(names(csv_data[i]))
  # not_ready_df <- csv_data[["stand_pep_quant_mergedcompomics.csv"]]
  # row_id <- grepl("peptide", colnames(csv_data[["stand_pep_quant_mergedcompomics.csv"]]))
  # cols_ids <- grepl("abundance", colnames(csv_data[["stand_pep_quant_mergedcompomics.csv"]]))
  # not_ready_df <- not_ready_df[c(1:No_rows),]
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

  new_df <- as.matrix(new_df[id_good, ])
  rownames(new_df) <- row_names
  return(new_df)
}

sizePortion <- function(input){
  cat("No of Rows: ", nrow(input), ",No of Cols: ",ncol(input), "\n", "nrow/ncol portion: ",
      round(digits = 2, nrow(input)/ncol(input)*100),"%")
  # return(nrow(input)/ncol(input))
}

# ### ** Examples
#
# # load part of the beatAML data
# beatAML_data <- NIMAA::beatAML[1:10000,]
#
# # convert to incidence matrix
# beatAML_incidence_matrix <- nominalAsBinet(beatAML_data)
#
# # extract submatrices with non-missing values
# sub_matrices <- extractSubMatrix(beatAML_incidence_matrix, col.vars = "patient_id",
#                                  row.vars = "inhibitor")



# read4_nimaa <- df_4_NIMAA(not_ready_df,"peptide", c("abundance_A", "abundance_B"))
read4_nimaa <- df_4_NIMAA(tmptbl, "protein_group", T, c("peptides_A", "peptides_B"))
cat("------------->>>>>>>>>>>>.  ",class(read4_nimaa))
subMat <- NULL
# extract submatrices with non-missing values
cat("--------> dataset name: ", nm, "-------> iteration: ", i, "\n")
matrixxx <- (read4_nimaa[1:1000, ])
sizePortion(t(matrixxx))
# rownames(matrixxx) <- NULL
# colnames(matrixxx) <- NULL
subMat <- extractSubMatrix(matrixxx, col.vars = "samples", row.vars = "protein")

extractSubMatrix(matrixxx, col.vars = NULL, row.vars = NULL)

rownames(matrixxx) <- rownames(read4_nimaa)[1:100]

# Call the function without meaningless labels
extractSubMatrix(matrixxx, col.vars = NULL, row.vars = NULL)

# x_after_arrange <- NIMAA::threeStepArrange(matrixxx)
# findSquareCutoffSubmatrix(x_after_arrange)



# load part of the beatAML data
# eg here
beatAML_data <- NIMAA::beatAML[1:794,]
# No of Rows:  10 ,No of Cols:  87
# nrow/ncol portion:  11.49 %

#  No of Rows:  11 ,No of Cols:  87
#  nrow/ncol portion:  12.64 %

# convert to incidence matrix
beatAML_incidence_matrix <- nominalAsBinet(beatAML_data)

sizePortion(beatAML_incidence_matrix)

# extract submatrices with non-missing values
# in following line the error is goning to be generated for the above size of the input dataset
sub_matrices <- extractSubMatrix(beatAML_incidence_matrix, col.vars = "patient_id",
                                 row.vars = "inhibitor")

# this is the error
# binmatnest.temperature
# 3.614421
# Error in `[.data.frame`(x, 1:i, 1:i) : undefined columns selected


beatAML_data <- NIMAA::beatAML[1:795,]

beatAML_data <- NIMAA::beatAML[1:10000,]
