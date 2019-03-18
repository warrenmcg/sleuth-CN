#' Find IQLR Denominators
#'
#' This function identifies which features in a matrix 
#' are suitable set to use as a denominator for the
#' 'Interquartile Logratio' (IQLR) transformation.
#' These are the features that have variance in the
#' interquartile range after transforming the data using
#' the centered logratio (CLR) transformation.
#' This was originally introduced in ALDEx2 in its paper.
#' 
#' @param mat a D x M matrix, with D features and M samples
#' @param base either "e" (for natural log scale) or "2" for base-2 logarithm
#'
#' @return a numeric vector containing the indices of the features that
#'   form the set to be used as the denominator for the IQLR transformation
#'
#' @references For the ALDEx2 paper discusing the IQLR transformation, see
#'   \url{https://dx.doi.org/10.1186\%2F2049-2618-2-15}
#' @importFrom stats quantile
#' @export
find_iqlr_denoms <- function(mat, base = "e") {
  clr_vals <- calculate_clr(mat, base)
  clr_var <- matrixStats::rowVars(clr_vals)
  qts <- stats::quantile(clr_var, na.rm = T)
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
#' @param remove_zeros boolean to see if this function
#'   should remove essential zeros (features with zeros in
#'   all samples). The default is \code{FALSE} to be
#'   compatible with sleuth, as its default filter removes
#'   essential zeros.
#' @param denom_method either 'geomean' or 'DESeq2' to
#'   use either the geometric mean of the IQLR features as the
#'   denominator, or the DESeq2-style size factors (focused on the median
#'   among the IQLR features) as the denominator. The IQLR features
#'   are selected using \code{\link{find_iqlr_denoms}}.
#' @param impute_method which method to use for imputing zeros.
#'   'multiplicative' (default) sets all values smaller than
#'   a imputation value 'delta' (determined by delta or
#'   impute_proportion) to that imputation value, and reduces
#'   all other values by the amount X * (1 - delta*num_zero_values /
#'   sum_constraint). 'additive' is similar to most other tools, and
#'   just adds the imputation value to all entries ('delta' must
#'   be specified)
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
iqlr_transformation <- function(mat, base = "e", remove_zeros = FALSE,
                                denom_method = "geomean",
                                impute_method = "multiplicative",
                                delta = NULL, impute_proportion = 0.65) {
  flip <- FALSE
  if (ncol(mat) > nrow(mat)) {
    mat <- t(mat)
    flip <- TRUE
  }

  if(remove_zeros) {
    mat <- remove_essential_zeros(mat)
  }

  imputed_mat <- impute_zeros(mat, method = impute_method, delta = delta,
                              impute_proportion = impute_proportion)
  denom_index <- find_iqlr_denoms(imputed_mat, base)
  iqlr_table <- calculate_alr(imputed_mat, base, denom_index, denom_method)
  if (flip) iqlr_table <- t(iqlr_table)
  iqlr_table
}
