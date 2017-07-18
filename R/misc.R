norm_identity <- function(mat) {
  rep(1, ncol(mat))
}

geomean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
