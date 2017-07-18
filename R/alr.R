#' Calculate ALR
#' 
#' This function calculates the additive logratio
#' transformation.
#'
#' @param mat an D x M matrix of D target IDs and
#'   M samples
#' @param denom_index a character vector of target IDs or
#'   an integer vector of row numbers. If there is
#'   more than one, the geometric mean of the 
#'   denominator values within one sample is used
#'   as the denominator.
#' 
#' @return (D - n) x M matrix of ALR-transformed
#'   values, with n equal to the number of denominator values 
calculate_alr <- function(mat, denom_index = NULL) {
  if(is.character(denom_index) & !(all(denom_index %in% rownames(mat))))
    stop(denom_index, " is not one of the row names of your matrix.")
  else if (!(all(denom_index %in% c(1:nrow(mat)))))
    stop(denom_index, " is outside of the bounds of your matrix.")

  if (any(mat == 0)) {
    stop("The ALR transformation cannot be done because there is ",
         "at least one zero value in the supplied matrix.")
  }

  if (length(denom_index) > 1) {
    denom_values <- mat[denom_index, ]
    denom_row <- apply(denom_values, 2, geomean)
  } else denom_row <- mat[denom_index, ]
  alr_table <- sweep(mat, 2, denom_row, "/")
  alr_table <- log(alr_table)
  alr_table <- alr_table[-denom_index, ]
  alr_table
}

#' Calculate ALR
#'
#' This function calculates the additive logratio
#' transformation. This function uses calculate_alr
#' internally after processing the data to handle
#' zeros. Essential zeros (i.e. target IDs with 
#' zeros in all samples) are excluded, and
#' round zeros (i.e. target IDs with at least one
#' zero value) are imputed using the mulitplicative
#' method.
#'
#' @param mat an D x M matrix of D target IDs and
#'   M samples
#' @param denom_name a character vector of target IDs or
#'   an integer vector of row numbers. If there is
#'   more than one, the geometric mean of the
#'   denominator values within one sample is used
#'   as the denominator.
#'
#' @return (D - n - z) x M matrix of ALR-transformed
#'   values, with n equal to the number of denominator values
#'   and z are the number of rows with essential zeros.
#'
#' @export
alr_transformation <- function(mat, denom_name = NULL) {
  stopifnot(is.character(denom_name))
  flip <- FALSE
  if (ncol(mat) > nrow(mat)) {
    mat <- t(mat)
    flip <- TRUE
  }
  mat <- remove_essential_zeros(mat)
  imputed_mat <- impute_rounded_zeros(mat)
  denom_index <- which(rownames(imputed_mat) == denom_name)
#  denom_values <- imputed_mat[denom_index, ]
  alr_table <- calculate_alr(imputed_mat, denom_index)
  if (flip) alr_table <- t(alr_table)
  alr_table
}
