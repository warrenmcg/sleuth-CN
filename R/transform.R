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
                            method = "multiplicative",
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
  e$method <- method
  e$base <- base
  if (type == "alr") {
    e$denom <- denom_name
    e$fun <- function(matrix, sf = 1, denom_name = eval(e$denom)) {
      alr_transformation(matrix, denom_name = denom_name,
                         method = e$method,
                         delta = e$delta,
                         impute_proportion = e$impute,
                         base = e$base)
    }
  } else if (type == "iqlr") {
    e$denom <- "iqlr"
    e$fun <- function(matrix, sf = 1) {
      iqlr_transformation(matrix,
                         delta = e$delta,
                         impute_proportion = e$impute,
                         base = e$base)
    }
  } else {
    e$denom <- "all"
    e$fun <- function(matrix, sf = 1) {
      clr_transformation(matrix,
                         delta = e$delta,
                         impute_proportion = e$impute,
                         base = e$base)
    }
  }
  transform_fun <- function(matrix, sf) { e$fun(matrix, sf) }
  return(transform_fun)
}

#' @export
get_norm_and_transform_funs <- function(type = "alr", denom_name = NULL,
                                        impute_method = "multiplicative",
                                        denom_method = "geomean",
                                        delta = NULL, impute_proportion = 0.65,
                                        base = "e") {
  base <- as.character(base)
  base <- match.arg(base, c("e", "2"))
  denom_method <- match.arg(denom_method, c("geomean", "DESeq2"))
  impute_method <- match.arg(impute_method, c("multiplicative", "additive"))
  type <- match.arg(type, c("alr", "ALR", "clr", "CLR", "iqlr", "IQLR"))
  type <- tolower(type)

  if (type == "alr" & is.null(denom_name)) {
    stop("you selected the 'ALR' transformation, but no ",
         "denominator was supplied. 'denom_name' is required ",
         "for the 'ALR' type.")
  }

  n <- new.env()
  if (type != "alr") {
    denom_name <- NULL
  }

  n$denoms <- denom_name
  n$type <- type
  n$method <- denom_method
  n$fun <- function(mat, denoms = n$denoms, type = n$type, method = n$method) {
    if (type == "clr") {
      denoms <- 1:nrow(mat)
    } else if (type == "iqlr") {
      denoms <- sleuthALR::find_iqlr_denoms(mat)
    } else if (type == "alr") {
      if (is(denoms, "character") && any(!(denoms %in% rownames(mat)))) {
        bad_denoms <- denoms[!(denoms %in% rownames(mat))]
        stop(paste("At least one of the supplied denominator features is",
                   "not found. Here is the list of denominators not found:",
                   paste(bad_denoms, collapse = ", ")))
      } else if (!is(denoms, "character")) {
        stop(paste("Class", class(denoms), "is unsupported to identify",
                   "denominator features for normalization"))
      }
    }

    if (method == "geomean") {
      if (length(denoms) == 1) {
        sf <- mat[denoms, ]
      } else {
        mat_nz <- mat[denoms, , drop = FALSE]
        sf <- apply(mat_nz, 2, sleuthALR::geomean)
      }
    } else {
      mat_nz <- mat[denoms, , drop = FALSE]
      p <- ncol(mat_nz)
      geo_means <- exp(rowMeans(log(mat_nz), na.rm = TRUE))
      s <- sweep(mat_nz, 1, geo_means, `/`)
      sf <- matrixStats::colMedians(s, na.rm = TRUE)
      scaling <- exp((-1/p) * sum(log(sf)))
      sf <- sf * scaling
    }
    sf
  }

  norm_func <- function(mat) n$fun(mat)
  environment(norm_func) <- n

  t <- new.env()
  t$delta <- delta
  t$impute <- impute_proportion
  t$method <- impute_method
  t$base <- base
  t$fun <- function(mat, sf = 1, delta = t$delta, impute = t$impute,
                    method = t$method, base = t$base) {
    logfunc <- switch(base, "e" = log, "2" = log2)
    mat <- suppressWarnings(
      sleuthALR:::impute_zeros(mat, delta = delta,
                               impute_proportion = impute,
                               method = method)
    )

    if (length(sf) == 1) {
      logfunc(mat / sf)
    } else if (length(sf) == ncol(mat)) {
      logfunc(t(t(mat) / sf))
    } else {
      stop("please provide a single size factor or a size factor vector equal ",
           "in length to the number of samples")
    }
  }
  transform_func <- function(mat, sf) t$fun(mat, sf = sf)
  environment(transform_func) <- t

  list(n_func = norm_func, t_func = transform_func)
}
