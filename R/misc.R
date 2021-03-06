#' Norm Identity
#'
#' Function to return size factors of just 1
#' for sleuth normalization (i.e. lack of
#' normalization). Normalization is not
#' necessary for logratio transformations.
#'
#' @param mat a numeric matrix, with columns for each
#'   sample, and rows for each feature
#' @return a vector of ones with length equal to the
#'   number of columns (samples) in the input matrix
#' @export
norm_identity <- function(mat) {
  rep(1, ncol(mat))
}

#' Geometric mean
#'
#' Function to return the geometric mean. By default,
#' it only calculates on positive numbers, and gives a
#' warning if there are any NA values or any numbers
#' that are zero or negative.
#'
#' @param x numeric vector of positive numbers.
#'
#' @return the geometric mean of x
#' @export
geomean <- function(x){
  if (any(is.na(x))) {
    num_na <- sum(is.na(x))
    warning(paste(num_na, "NA values were found. Omitting them."))
  }
  if (length(which(x <= 0)) > 0) {
    num_zero <- length(which(x <= 0))
    warning(paste(num_zero, "values were zero or negative. Omitting them."))
  }
  x_vals <- x[!is.na(x) & x > 0]
  exp(sum(log(x_vals)) / length(x_vals))
}

#' Get DESeq2 Size Factors
#'
#' Function to return the DESeq2-style size factors of a collection of features.
#' The geometric mean is calculated for feature across samples. Then the
#' ratio of every sample observation to the geometric mean across samples is
#' calculated. From there, the median of ratios within each sample is selected
#' as the size factor.
#'
#' @param x a D x M matrix with D features (the denominator features) and
#'   M samples
#' @param denoms a vector of either indices or feature names that match row names
#'   of the matrix. These are the ones that will be used for calculating the
#'   compositional normalization using the DESeq2-method.
#'
#' @return a vector of the size factors for each sample
#' @references For a discussion of why the DESeq2 size factor is compatible
#'   with Compositional Normalization, see
#'   \url{https://dx.doi.org/10.1093/bioinformatics/bty175}
#' @references For the original DESeq2 method, see
#'   \link[DESeq2]{estimateSizeFactorsForMatrix}
#' @importFrom methods is
#' @importFrom stats median
#' @export
deseq_size_factors <- function(x, denoms = NULL) {
  stopifnot(is(x, "matrix"))
  if (!is.null(denoms)) {
    if (is(denoms, "character") && !all(denoms %in% rownames(x))) {
      missing_rows <- denoms[which(!denoms %in% rownames(x))]
      formatted_rows <- paste(missing_rows, collapse = "', '")
      err_msg <- paste0("There were ", length(missing_rows), " features ",
                        "from the rownames of the supplied matrix. Here are ",
                        "missing features: '", formatted_rows, "'")
      stop(err_msg)
    }
    x <- x[denoms, , drop = FALSE]
  }
  logmeans <- rowMeans(log(x))
  sf <- apply(x, 2, function(sample) {
    ratios <- (log(sample) - logmeans)[is.finite(logmeans) & sample > 0]
    exp(median(ratios))
  })
  sf
}

#' Get Denom Name(s)
#' 
#' Function to retrieve the denominator(s) for a sleuth object
#' processed using compositional data analysis. If the analysis
#' was done using CLR (centered logratio), the denominator name is
#' "all".
#'
#' @param obj the sleuth object
#' @return The name of the denominator(s) for the sleuth object
#' @export
get_denom_names <- function(obj) {
  stopifnot(methods::is(obj, "sleuth"))
  name <- environment(obj$norm_fun_counts)$denom_name
  if (is.null(name)) {
    name <- environment(obj$transform_fun_counts)$denom_name
    if (is.null(name)) {
      name <- environment(obj$transform_fun)$denom_name
      if (is.null(name)) {
        stop("This sleuth object was not processed using compositional data analysis.",
             "\nThus its denominator does not exist.")
      }
    }
  }

  name
}

# This function retrieves the proper weighting function for 'sleuth_alr_results'.
# It is exponentiating the mean observations using the inverse of the logarithm
# used for the original transformation
get_alr_weight <- function(obj) {
  stopifnot(methods::is(obj, "sleuth"))
  base <- environment(obj$transform_fun)$base
  base <- as.character(base)
  if (base == "e")
    func <- exp
  else (base == "2")
    func <- function(x) 2^x
  func
}

# function to clean the sleuth object to remove the denominator names
# from the processed data to prevent the denominator from affecting the modeling
# steps
clean_denom_names <- function(obj, lr_method) {
  stopifnot(methods::is(obj, 'sleuth'))
  denom_names <- get_denom_names(obj)
  len <- length(denom_names)
  if (len > 1 || (is.null(denom_names) || denom_names %in% c("all", "iqlr"))) {
    return(obj)
  } else {
    if(!is.null(obj$bs_summary$obs_counts)) {
      indices <- which(rownames(obj$bs_summary$obs_counts) %in% denom_names)
      obj$bs_summary$obs_counts <- obj$bs_summary$obs_counts[-indices, ]
      obj$bs_summary$sigma_q_sq <- obj$bs_summary$sigma_q_sq[-indices]
    }
    if(!is.null(obj$bs_summary$obs_tpm)) {
      indices <- which(rownames(obj$bs_summary$obs_tpm) %in% denom_names)
      obj$bs_summary$obs_tpm <- obj$bs_summary$obs_tpm[-indices, ]
      obj$bs_summary$sigma_q_sq_tpm <- obj$bs_summary$sigma_q_sq_tpm[-indices]
    }
    # If both norm and transform functions are used, then the denominator
    # feature will still have inferential variation, so its bs_quants can
    # be used. For 'transform' only, though, each bootstrap is normalized
    # so the inferential variation will also be zero.
    if (lr_method == "transform") {
      if(length(obj$bs_quants) > 0 & length(obj$bs_quants[[1]]) > 0) {
        obj$bs_quants <- lapply(obj$bs_quants, function(sample) {
          lapply(sample, function(summary) {
            indices <- which(rownames(summary) %in% denom_names)
            summary[-indices, ]
          })
        })
      }
    }
    return(obj)
  }
}

# function copied from sleuth to get a test. This is used for 'sleuth_alr_results'.
get_alr_test <- function (obj, label, type, model)
{
    stopifnot(methods::is(obj, "sleuth"))
    stopifnot(type %in% c("lrt", "wt"))
    res <- NULL
    if (type == "lrt") {
        res <- obj$tests[[type]][[label]]
    }
    else {
        if (missing(model)) {
            stop("must specify a model with wald test")
        }
        res <- obj$tests[[type]][[model]][[label]]
    }
    if (is.null(res)) {
        stop("'", label, "' is not a valid label for a test.",
            " Please see valid models and tests using the functions 'models' and 'tests'.",
            " Remember to also correctly specify the test type.")
    }
    res
}
