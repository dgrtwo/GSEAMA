#' Standard-evaluation version of sparse_cast
#'
#' @param data A tbl
#' @param row_col String version of column to use as row names
#' @param column_col String version of column to use as column names
#' @param value_col String version of column to use as sparse matrix values,
#' or a numeric value to use
#'
#' @import Matrix
sparse_cast_ <- function(data, row_col, column_col, value_col = NULL) {
  row_names <- data[[row_col]]
  col_names <- data[[column_col]]
  if (is.numeric(value_col)) {
    values <- value_col
  } else {
    values <- data[[value_col]]
  }
  
  # if it's a factor, preserve ordering
  row_u <- if (is.factor(row_names)) levels(row_names) else unique(row_names)
  col_u <- if (is.factor(col_names)) levels(col_names) else unique(col_names)
  
  ret <- Matrix(0, nrow = length(row_u), ncol = length(col_u),
                dimnames = list(as.character(row_u), as.character(col_u)),
                sparse = TRUE)
  
  ret[cbind(match(row_names, row_u), match(col_names, col_u))] <- values
  
  ret
}


#' Create a sparse matrix from row names, column names, and values
#' in a table.
#'
#' @param data A tbl
#' @param row A bare column name to use as row names in sparse matrix
#' @param column A bare column name to use as column names in sparse matrix
#' @param value A bare column name to use as sparse matrix values, default 1
#'
#' @return A sparse Matrix object, with one row for each unique value in
#' the \code{row} column, one column for each unique value in the \code{column}
#' column, and with as many non-zero values as there are rows in \code{data}.
#'
#' @examples
#'
#' dat <- data.frame(a = c("row1", "row1", "row2", "row2", "row2"),
#'   b = c("col1", "col2", "col1", "col3", "col4"),
#'   val = 2)
#'
#' sparse_cast(dat, a, b)
#'
#' sparse_cast(dat, a, b, val)
#'
#' @name sparse_cast
sparse_cast <- function(data, row, column, value) {
  if (missing(value)) {
    value_col <- 1
  } else {
    value_col <- as.character(substitute(value))
    if (is.null(data[[value_col]])) {
      value_col <- value
    }
  }
  
  sparse_cast_(data, as.character(substitute(row)),
               as.character(substitute(column)), value_col)
}


#' Trim the length of a character vector, adding ellipses to trimmed elements
#' 
#' This utility function is used to shorten column names so they can appear
#' on a figure.
#' 
#' @param x character vector
#' @param width length of string to trim to
#' 
#' @export
trim_ellipses <- function(x, width) {
  ret <- stringr::str_sub(x, 1, width)
  toolong <- stringr::str_length(x) > width
  ret[toolong] <- stringr::str_c(ret[toolong], "...")
  ret
}
