#' Get Logration Function
#'
#' This helps you create a function that can be input as a
#' transformation function for 'sleuth_prep'.
#' 
#' @param type either "alr", "clr", or "iqlr" (also accepts "ALR" / "CLR" / "IQLR")
#' @param denom_name a character vector of target IDs that
#'  will be the denominator. More than one is acceptable, and
#'  you can also input a numeric vector of row numbers instead.
#'  Note that this will be ignored if 'clr' or 'iqlr' is used.
#' @param delta a number that is the imputed value. If \code{NULL},
#'  delta = impute_proportion * (minimum value in sample) 
#' @param impute_proportion percentage of minimum value that
#'  becomes the imputed value. Only used if delta is \code{NULL}
#' @param base the base used for the logarithm. Currently only supports
#'  "e" or "2" (can also specify the number 2).
#'
#' @return a transformation function ready to be used for the
#'  'transform_fun' option in 'sleuth_prep'
#' 
#' @export
get_lr_function <- function(type = "alr", denom_name = NULL,
                            delta = NULL, impute_proportion = 0.65,
                            base = "e") {
  type <- match.arg(type, c("alr", "ALR", "clr", "CLR", "iqlr", "IQLR"))
  type <- tolower(type)
  
  if (type == "alr" & is.null(denom_name)) {
    stop("you selected the 'ALR' transformation, but no ",
         "denominator was supplied. 'denom_name' is required ",
         "for the 'ALR' type.")
  }
  e <- new.env()
  e$delta <- delta
  e$impute <- impute_proportion
  e$base <- base
  if (type == "alr") {
    e$denom <- denom_name
    e$fun <- function(matrix, denom_name = eval(e$denom)) {
      alr_transformation(matrix, denom_name = denom_name,
                         delta = e$delta,
                         impute_proportion = e$impute,
                         base = e$base)
    }
  } else if (type == "iqlr") {
    e$denom <- "iqlr"
    e$fun <- function(matrix) {
      iqlr_transformation(matrix,
                         delta = e$delta,
                         impute_proportion = e$impute,
                         base = e$base)
    }
  } else {
    e$denom <- "all"
    e$fun <- function(matrix) {
      clr_transformation(matrix,
                         delta = e$delta,
                         impute_proportion = e$impute,
                         base = e$base)
    }
  }
  transform_fun <- function(matrix) { e$fun(matrix) }
  return(transform_fun)
}
