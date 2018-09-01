#' Impute Rounded Zeros
#' 
#' This uses various strategies to impute zeros
#' that are considered "rounded zeros", that is, values
#' that are zero because the true value is positive but below
#' the detection limit (i.e. the sequencing depth). See the
#' references for more information.
#'
#' @param mat a numeric matrix to be imputed, with N
#'   targets and M samples. This method assumes that N > M.
#'   N x M and M x N matrices are both valid inputs.
#' @param method the choice of how to impute the rounded zeros.
#'   only "multiplicative" and "additive" is supported at this time.
#' @param delta the value to impute; if NULL,
#'   delta = impute_proportion x (the minimum value in a sample)
#' @param impute_proportion the proportion of the minimum value
#'   that should be the delta value (default = 65\%).
#'
#' @return a matrix with the same dimensions as input and all zero
#'   values imputed by sample.
#' @export
impute_zeros <- function(mat, method = "multiplicative",
                         delta = NULL, impute_proportion = 0.65) {
  flip <- FALSE
  if (ncol(mat) > nrow(mat)) {
    mat <- t(mat)
    flip <- TRUE
  }

  method <- match.arg(method, c("multiplicative", "additive", "EM", "bayesian"))
  stopifnot(impute_proportion > 0 & impute_proportion < 1)
  stopifnot(is.null(delta) | delta > 0)

  if (any(rowSums(mat) == 0)) {
    bad_rows <- which(rowSums(mat) == 0)
    num_bad <- length(bad_rows)
    # Unless you are normalizing sleuth bootstraps, which may have zeros in all samples
    warning("There were ", num_bad, " targets that had all zero values. ",
            "Those targets should be removed using 'remove_essential_zeros'")
  }

  if (method == "multiplicative") {
    tmp_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    if(is.null(delta)) {
      bool <- mat >= .Machine$double.eps / impute_proportion
      tmp_mat[bool] <- mat[bool]
      mat[!bool] <- .Machine$double.eps / impute_proportion
      detection_limit <- suppressWarnings(exp(matrixStats::colMins(log(mat), na.rm = T)))
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
  } else if (method == "additive") {
    stopifnot(!is.null(delta) && delta > 0)
    new_mat <- mat + delta
  } else {
    stop("Methods 'EM' and 'bayesian' are not supported yet.")
  }

  if (flip) {
    new_mat <- t(new_mat)
  }
  new_mat
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
#' @export
remove_essential_zeros <- function(mat) {
  allZeroRows <- matrixStats::rowAlls(mat, value = 0)
  mat[-allZeroRows, ]
}
