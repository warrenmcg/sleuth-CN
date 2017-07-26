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
