find_iqlr_denoms <- function(mat, base = "e") {
  clr_vals <- calculate_clr(mat, base)
  clr_var <- apply(clr_vals, 1, var)
  qts <- quantile(clr_var, na.rm = T)
  indices <- which(clr_var > qts[2] & clr_var < qts[4])
  indices
}

#' Calculate IQLR
#'
#' This function calculates the interquartile logratio
#' transformation. This function first processes the data
#' to handle zeros. It then calculates which features
#' have the variance across samples within the interquartile range
#' after the CLR transformation. It then uses calculate_alr
#' internally after processing the data to handle
#' zeros. Essential zeros (i.e. target IDs with 
#' zeros in all samples) are excluded, and
#' round zeros (i.e. target IDs with at least one
#' zero value) are imputed using the mulitplicative
#' method.
#'
#' @param mat an D x M matrix of D target IDs and
#'   M samples
#' @param base what should the base of the logarithm be?
#'   currently only supports base "e" and base 2.
#' @param delta a number that is the imputed value. If \code{NULL},
#'  delta = impute_proportion * (minimum value in sample)
#' @param impute_proportion percentage of minimum value that
#'  becomes the imputed value. Only used if delta is \code{NULL}
#' 
#' @return (D - n - z) x M matrix of ALR-transformed
#'   values, with n equal to the number of denominator values
#'   and z are the number of rows with essential zeros.
#'
#' @export
iqlr_transformation <- function(mat, base = "e",
                                remove_zeros = FALSE, delta = NULL,
                                impute_proportion = 0.65) {
  flip <- FALSE
  if (ncol(mat) > nrow(mat)) {
    mat <- t(mat)
    flip <- TRUE
  }

  if(remove_zeros) {
    mat <- remove_essential_zeros(mat)
  }

  imputed_mat <- impute_zeros(mat, delta = delta,
                                      impute_proportion = impute_proportion)
  denom_index <- find_iqlr_denoms(imputed_mat, base)
  iqlr_table <- calculate_alr(imputed_mat, base, denom_index)
  if (flip) iqlr_table <- t(iqlr_table)
  iqlr_table
}
