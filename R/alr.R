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
#' @param denom_method either 'geomean' or 'DESeq2' to
#'   use either the geometric mean of the features as the
#'   denominator, or the DESeq2-style size factors
#'   as the denominator. If there is only one feature,
#'   'geomean' uses the feature's value itself, whereas
#'   'DESeq2' will use the ratio of the feature to 
#'   the geometric mean across samples.
#' 
#' @return D x M matrix of ALR-transformed values
calculate_alr <- function(mat, base = "e", denom_index = NULL, denom_method = "geomean") {
  if(is.character(denom_index) & !(all(denom_index %in% rownames(mat))))
    stop(denom_index, " is not one of the row names of your matrix.")
  else if (!(all(denom_index %in% c(1:nrow(mat)))))
    stop(denom_index, " is outside of the bounds of your matrix.")

  base <- as.character(base)
  base <- match.arg(base, c("e", "2"))
  denom_method <- match.arg(denom_method, c("geomean", "DESeq2"))

  if (any(mat == 0)) {
    stop("The ALR transformation cannot be done because there is ",
         "at least one zero value in the supplied matrix.")
  }

  if (denom_method == "geomean") {
    if (length(denom_index) > 1) {
      denom_values <- mat[denom_index, ]
      denom_row <- apply(denom_values, 2, geomean)
    } else denom_row <- mat[denom_index, ]
  } else {
    denom_row <- deseq_size_factors(mat, denoms = denom_index)
  }

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
#' zeros in all samples) can be optionally excluded, and
#' all remaining zeros, including rounded zeros
#' (i.e. target IDs with at least one zero value) are
#' always imputed using the hidden impute_zeros function.
#'
#' @param mat an D x M matrix of D target IDs and
#'   M samples
#' @param denom_name a character vector of target IDs or
#'   an integer vector of row numbers. If there is
#'   more than one, the geometric mean of the
#'   denominator values within one sample is used
#'   as the denominator.
#' @param base what should the base of the logarithm be?
#'   currently only supports base "e" and base 2.
#' @param remove_zeros boolean to see if this function
#'   should remove essential zeros (features with zeros in
#'   all samples). The default is \code{FALSE} to be
#'   compatible with sleuth, as its default filter removes
#'   essential zeros.
#' @param denom_method either 'geomean' or 'DESeq2' to
#'   use either the geometric mean of the features as the
#'   denominator, or the DESeq2-style size factors
#'   as the denominator. If there is only one feature,
#'   'geomean' uses the feature's value itself, whereas
#'   'DESeq2' will use the ratio of the feature to 
#'   the geometric mean across samples.
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
alr_transformation <- function(mat, denom_name,
                               base = "e", remove_zeros = FALSE,
                               denom_method = "geomean",
                               impute_method = "multiplicative",
                               delta = NULL,
                               impute_proportion = 0.65) {
  stopifnot(methods::is(mat, "matrix"))
  stopifnot(methods::is(mat[,1], "numeric"))

  flip <- FALSE
  if (ncol(mat) > nrow(mat)) {
    mat <- t(mat)
    flip <- TRUE
  }

  if (remove_zeros) {
    mat <- remove_essential_zeros(mat)
  }

  imputed_mat <- impute_zeros(mat, method = impute_method, delta = delta,
                              impute_proportion = impute_proportion)
  denom_index <- ifelse(is.character(denom_name),
    denom_index <- which(rownames(imputed_mat) %in% denom_name),
    denom_index <- denom_name
  )
  alr_table <- calculate_alr(imputed_mat, base, denom_index, denom_method)
  if (flip) alr_table <- t(alr_table)
  alr_table
}
