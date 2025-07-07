#  NIMAA inspection

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

beatAML_incidence_matrix <- nominalAsBinet(beatAML_data)

sizePortion(beatAML_incidence_matrix)

# extract submatrices with non-missing values
# in following line the error is goning to be generated for the above size of the input dataset
sub_matrices <- extractSubMatrix(beatAML_incidence_matrix, col.vars = "patient_id",
                                 row.vars = "inhibitor")
