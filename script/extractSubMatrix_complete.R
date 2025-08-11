ns <- asNamespace("NIMAA")  # replace with the correct package name

threeStepArrange            <- get("threeStepArrange", envir = ns)
findElementMaxSubmatrix     <- get("findElementMaxSubmatrix", envir = ns)
findSquareCutoffSubmatrix   <- get("findSquareCutoffSubmatrix", envir = ns)
findRowwiseCutoffSubmatrix  <- get("findRowwiseCutoffSubmatrix", envir = ns)
findColwiseCutoffSubmatrix  <- get("findColwiseCutoffSubmatrix", envir = ns)
plotSubmatrix               <- get("plotSubmatrix", envir = ns)

extractSubMatrix_complete <- function (x, shape = "All", verbose = FALSE, palette = "Greys",
          row.vars = NULL, col.vars = NULL, bar = 1, plot_weight = FALSE,
          print_skim = FALSE)
{
  input_type <- c("matrix", "data.frame")
  if (!(class(x)[1] %in% input_type)) {
    stop("The input should be matrix or data frame")
  }
  req_shape <- c("Square", "Rectangular_row", "Rectangular_col",
                 "Rectangular_element_max", "All")
  if (FALSE %in% (shape %in% req_shape)) {
    stop("Please specify the desire shape of the output")
  }
  R.names <- rownames(x)
  C.names <- colnames(x)
  if (is.null(R.names) | is.null(C.names)) {
    warning("There is no row names or column names for given input")
  }
  if (is.null(row.vars) | is.null(col.vars)) {
    warning("There is no name specified for samples on x-axis lable or y-axis lable")
  }
  if (is.null(row.vars) | is.null(col.vars)) {
    c_lable <- "columns"
    r_lable <- "rows"
  }
  else {
    c_lable <- as.character(col.vars)
    r_lable <- as.character(row.vars)
  }
  x_with0 <- x
  x_with0[is.na(x_with0)] <- 0
  binmatnest2.temp <- bipartite::nested(x_with0)
  print(binmatnest2.temp)
  if (binmatnest2.temp < 1) {
    cat(crayon::red$bold("The nestedness temperature is less than 1, highly nested!\nWe suggest that divide the data into different parts."))
    if (!utils::askYesNo("Do you still want to continue?")) {
      stop("Stop the funtion by user.")
    }
  }
  x <- x %>% as.data.frame()
  t_x <- t(x) %>% as.data.frame()
  x_after_arrange_without_weight <- threeStepArrange(x)
  x_after_arrange_with_weight <- x[rownames(x_after_arrange_without_weight),
                                   colnames(x_after_arrange_without_weight)]
  t_x_after_arrange_without_weight <- threeStepArrange(t_x)
  t_x_after_arrange_with_weight <- t_x[rownames(t_x_after_arrange_without_weight),
                                       colnames(t_x_after_arrange_without_weight)]
  result <- list()
  rect_data <- list()
  rect_data_t <- list()
  if ("Square" %in% shape || "All" %in% shape) {
    max_sq_nomiss <- findSquareCutoffSubmatrix(x_after_arrange_without_weight,
                                               bar)
    max_sq_nomiss_t <- findSquareCutoffSubmatrix(t_x_after_arrange_without_weight,
                                                 bar)
    rect_data$max_sq_nomiss <- max_sq_nomiss
    rect_data_t$max_sq_nomiss_t <- max_sq_nomiss_t
    x[rownames(max_sq_nomiss), colnames(max_sq_nomiss), drop = FALSE]
    if (prod(dim(max_sq_nomiss)) < prod(dim(max_sq_nomiss_t))) {
      max_sq_nomiss_better <- t(max_sq_nomiss_t)
    }
    else {
      max_sq_nomiss_better <- max_sq_nomiss
    }
    cat(crayon::green$bold("Size of Square: \t", dim(max_sq_nomiss_better)[1],
                           "rows x ", dim(max_sq_nomiss_better)[2], "columns",
                           "\n"))
    max_sq_nomiss_better <- x[rownames(max_sq_nomiss_better),
                              colnames(max_sq_nomiss_better), drop = FALSE]
    result$Square <- max_sq_nomiss_better
  }
  if ("Rectangular_row" %in% shape || "All" %in% shape) {
    row_max_rec <- findRowwiseCutoffSubmatrix(x_after_arrange_without_weight,
                                              bar)
    row_max_rec_t <- findColwiseCutoffSubmatrix(t_x_after_arrange_without_weight,
                                                bar)
    rect_data$row_max_rec <- row_max_rec
    rect_data_t$row_max_rec_t <- row_max_rec_t
    if (prod(dim(row_max_rec)) < prod(dim(row_max_rec_t))) {
      row_max_rec_better <- t(row_max_rec_t)
    }
    else {
      row_max_rec_better <- row_max_rec
    }
    cat(crayon::green$bold("Size of Rectangular_row: \t",
                           dim(row_max_rec_better)[1], "rows x ", dim(row_max_rec_better)[2],
                           "columns", "\n"))
    row_max_rec_better <- x[rownames(row_max_rec_better),
                            colnames(row_max_rec_better), drop = FALSE]
    result$Rectangular_row <- row_max_rec_better
  }
  if ("Rectangular_col" %in% shape || "All" %in% shape) {
    col_max_rec <- findColwiseCutoffSubmatrix(x_after_arrange_without_weight,
                                              bar)
    col_max_rec_t <- findRowwiseCutoffSubmatrix(t_x_after_arrange_without_weight,
                                                bar)
    rect_data$col_max_rec <- col_max_rec
    rect_data_t$col_max_rec_t <- col_max_rec_t
    if (prod(dim(col_max_rec)) < prod(dim(col_max_rec_t))) {
      col_max_rec_better <- t(col_max_rec_t)
    }
    else {
      col_max_rec_better <- col_max_rec
    }
    cat(crayon::green$bold("Size of Rectangular_col: \t",
                           dim(col_max_rec_better)[1], "rows x ", dim(col_max_rec_better)[2],
                           "columns", "\n"))
    col_max_rec_better <- x[rownames(col_max_rec_better),
                            colnames(col_max_rec_better), drop = FALSE]
    result$Rectangular_col <- col_max_rec_better
  }
  if ("Rectangular_element_max" %in% shape || "All" %in% shape) {
    max_rec <- findElementMaxSubmatrix(x_after_arrange_without_weight)
    max_rec_t <- findElementMaxSubmatrix(t_x_after_arrange_without_weight)
    rect_data$max_rec <- max_rec
    rect_data_t$max_rec <- max_rec_t
    if (prod(dim(max_rec)) < prod(dim(max_rec_t))) {
      max_rec_better <- t(max_rec_t)
    } else {
      max_rec_better <- max_rec
    }
    cat(crayon::green$bold(
      "Size of Rectangular_element_max: \t",
      dim(max_rec_better)[1], "rows x ", dim(max_rec_better)[2],
      "columns", "\n"
    ))
    max_rec_better <- x[rownames(max_rec_better), colnames(max_rec_better), drop = FALSE]

    ## --- NEW: Find dropped parts ---
    dropped_rows <- setdiff(rownames(x_after_arrange_without_weight), rownames(max_rec_better))
    dropped_cols <- setdiff(colnames(x_after_arrange_without_weight), colnames(max_rec_better))
    dropped_part <- x_after_arrange_without_weight[dropped_rows, dropped_cols, drop = FALSE]

    ## Store both kept and dropped parts
    result$Rectangular_element_max <- max_rec_better
    result$Rectangular_element_max_dropped <- dropped_part
  }
  if (plot_weight) {
    x_after_arrange <- x_after_arrange_with_weight
    t_x_after_arrange <- t_x_after_arrange_with_weight
  }
  else {
    x_after_arrange <- x_after_arrange_without_weight
    t_x_after_arrange <- t_x_after_arrange_without_weight
  }
  print_dataframe <- x_after_arrange %>% as.matrix() %>% as.data.frame.table()
  print_dataframe_t <- t_x_after_arrange %>% as.matrix() %>%
    as.data.frame.table()
  plotSubmatrix(x = x, print_dataframe = print_dataframe, rect_data = rect_data,
                verbose = verbose, palette = palette, c_lable = c_lable,
                r_lable = r_lable, figure_name = "Row_wise_arrangement",
                title = "Row-wise-arrangement")
  plotSubmatrix(x = t_x, print_dataframe = print_dataframe_t,
                rect_data = rect_data_t, verbose = verbose, palette = palette,
                c_lable = r_lable, r_lable = c_lable, figure_name = "Column_wise_arrangement",
                title = "Column-wise-arrangement")
  if (print_skim) {
    for (mat in names(result)) {
      switch(mat, Rectangular_row = {
        Rectangular_row <- result$Rectangular_row
        print(skimr::skim(Rectangular_row))
      }, Rectangular_element_max = {
        Rectangular_element_max <- result$Rectangular_element_max
        print(skimr::skim(Rectangular_element_max))
      }, Square = {
        Square <- result$Square
        print(skimr::skim(Square))
      }, Rectangular_col = {
        Rectangular_col <- result$Rectangular_col
        print(skimr::skim(Rectangular_col))
      })
    }
  }
  return(result)
}
