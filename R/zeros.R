#' Impute Rounded Zeros
#' 
#' This uses the multiplicative strategy to impute zeros
#' that are considered "rounded zeros", that is, values
#' that are zero because the true value is positive but below
#' the detection limit (i.e. the sequencing depth). See the
#' references for more information.
#'
#' @param mat an N x M numeric matrix to be imputed, with N
#'   targets and M samples.
#' @param method the choice of how to impute the rounded zeros.
#'   only "multiplicative" is supported at this time.
#' @param delta the value to impute; if NULL,
#'   delta = impute_proportion x (the minimum value in a sample)
#' @param impute_proportion the proportion of the minimum value
#'   that should be the delta value (default = 65\%).
#'
#' @return an N x M matrix with all zero values imputed.
#'
impute_zeros <- function(mat, method = "multiplicative",
                         delta = NULL, impute_proportion = 0.65) {
  method <- match.arg(method, c("multiplicative", "EM", "bayesian"))
  stopifnot(impute_proportion > 0 & impute_proportion < 1)
  stopifnot(is.null(delta) | delta > 0)
  if (method == "multiplicative") {
    tmp_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    if(is.null(delta)) {
      bool <- mat >= .Machine$double.eps / impute_proportion
      tmp_mat[bool] <- mat[bool]
      detection_limit <- suppressWarnings(exp(matrixStats::colMins(log(tmp_mat), na.rm = T)))
      delta <- detection_limit * impute_proportion
    } else {
      bool <- mat >= delta
      tmp_mat[bool] <- mat[bool]
    }
    c_i <- colSums(tmp_mat)
    zeros <- tmp_mat == 0
    num_zeros <- colSums(zeros)
    new_mat <- sweep(tmp_mat, 2, (1 - (delta * num_zeros) / c_i), "*")
    new_mat[zeros] <- delta
    new_mat
  } else {
    stop("Methods 'EM' and 'bayesian' are not supported yet.")
  }
}

#' Remove Essential Zeros
#' 
#' This function excludes "essential zeros", that is,
#' target IDs that are assumed to be silent and
#' non-expressed.
#' 
#' @param mat an N x M numeric matrix.
#' 
#' @return an (N-z) x M numeric matrix, with z
#'   equal to the number of rows with all zero values.
remove_essential_zeros <- function(mat) {
  allZeroRows <- matrixStats::rowAlls(mat, value = 0)
  mat[-allZeroRows, ]
}
