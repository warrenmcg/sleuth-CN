#' Calculate ALR
#' 
#' This function calculates the additive logratio
#' transformation.
#'
#' @param mat an D x M matrix of D target IDs and
#'   M samples
#' @param base what should the base of the logarithm be?
#'   currently only supports base "e" and base 2.
#' @param denom_index a character vector of target IDs or
#'   an integer vector of row numbers. If there is
#'   more than one, the geometric mean of the 
#'   denominator values within one sample is used
#'   as the denominator.
#' 
#' @return (D - n) x M matrix of ALR-transformed
#'   values, with n equal to the number of denominator values 
calculate_alr <- function(mat, base = "e", denom_index = NULL) {
  if(is.character(denom_index) & !(all(denom_index %in% rownames(mat))))
    stop(denom_index, " is not one of the row names of your matrix.")
  else if (!(all(denom_index %in% c(1:nrow(mat)))))
    stop(denom_index, " is outside of the bounds of your matrix.")

  base <- as.character(base)
  base <- match.arg(base, c("e", "2"))

  if (any(mat == 0)) {
    stop("The ALR transformation cannot be done because there is ",
         "at least one zero value in the supplied matrix.")
  }

  if (length(denom_index) > 1) {
    denom_values <- mat[denom_index, ]
    denom_row <- apply(denom_values, 2, geomean)
  } else denom_row <- mat[denom_index, ]
  alr_table <- sweep(mat, 2, denom_row, "/")
  if (base == "e") alr_table <- log(alr_table) else
    alr_table <- log(alr_table, as.integer(base))
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
#' @param base what should the base of the logarithm be?
#'   currently only supports base "e" and base 2.
#' @param denom_name a character vector of target IDs or
#'   an integer vector of row numbers. If there is
#'   more than one, the geometric mean of the
#'   denominator values within one sample is used
#'   as the denominator.
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
alr_transformation <- function(mat, base = "e",
                               denom_name = NULL, delta = NULL,
                               impute_proportion = 0.65) {
#  stopifnot(is.character(denom_name))
  flip <- FALSE
  if (ncol(mat) > nrow(mat)) {
    mat <- t(mat)
    flip <- TRUE
  }
#  mat <- remove_essential_zeros(mat)
  imputed_mat <- impute_rounded_zeros(mat, delta = delta,
                                      impute_proportion = impute_proportion)
  if (is.character(denom_name))
    denom_index <- which(rownames(imputed_mat) == denom_name)
  else
    denom_index <- denom_name
  alr_table <- calculate_alr(imputed_mat, base, denom_index)
  if (flip) alr_table <- t(alr_table)
  alr_table
}
