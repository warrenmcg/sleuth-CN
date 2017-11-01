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

#' Get Denom Name
#' 
#' Function to retrieve the denominator for a sleuth object
#' processed using compositional data analysis.
#'
#' @param obj the sleuth object
#' @return The name of the denominator(s) for the sleuth object
#' @export
get_denom_name <- function(obj) {
  stopifnot(is(obj, "sleuth"))
  name <- environment(obj$transform_fun)$denom_name
  if (is.null(name)) {
    stop("This sleuth object was not processed using compositional data analysis.",
         "\nThus its denominator does not exist.")
  }
  name
}
