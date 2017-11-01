#' Norm Identity
#'
#' Function to return size factors of just 1
#' for sleuth normalization (i.e. lack of
#' normalization). Normalization is not
#' necessary for logratio transformations.
#'
#' @export
norm_identity <- function(mat) {
  rep(1, ncol(mat))
}

geomean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

#' Get Denom Name(s)
#' 
#' Function to retrieve the denominator(s) for a sleuth object
#' processed using compositional data analysis. If the analysis
#' was done using CLR (centered logratio), the denominator name is
#' all.
#'
#' @param obj the sleuth object
#' @return The name of the denominator(s) for the sleuth object
#' @export
get_denom_names <- function(obj) {
  stopifnot(is(obj, "sleuth"))
  name <- environment(obj$transform_fun)$denom_name
  if (is.null(name)) {
    stop("This sleuth object was not processed using compositional data analysis.",
         "\nThus its denominator does not exist.")
  }
  name
}

# function to clean the sleuth object to remove the denominator names
# from the processed data to prevent the denominator from affecting the modeling
# steps
clean_denom_names <- function(obj) {
  stopifnot(is(obj, 'sleuth'))
  denom_names <- get_denom_names(obj)
  if (is.null(denom_names) | denom_names == "all") {
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
    if(length(obj$bs_quants) > 0 & length(obj$bs_quants[[1]]) > 0) {
      obj$bs_quants <- lapply(obj$bs_quants, function(sample) {
        lapply(sample, function(summary) {
          indices <- which(rownames(summary) %in% denom_names)
          summary[-indices, ]
        })
      })
    }
    return(obj)
  }
}
