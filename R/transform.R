#' Get Logratio Function
#'
#' This helps you create a function that can be input as a
#' transformation function for 'sleuth_prep'.
#' 
#' @param type either "alr", "clr", or "iqlr" (also accepts "ALR" / "CLR" / "IQLR")
#' @param denom_name a character vector of target IDs that
#'  will be the denominator. More than one is acceptable, and
#'  you can also input a numeric vector of row numbers instead.
#'  Note that this will be ignored if 'clr' or 'iqlr' is used.
#' @param lr_method the choice of how to conduct compositional normalization.
#'   "both" provides a compositional normalization function and a compositional
#'   transformation functions;
#'   "transform" provides a compositional transformation function that also does the normalization
#' @param denom_method the choice of what kind of compositional normalization to do
#'   when more than one feature is used. "geomean" takes the geometric of all features within a sample
#'   as the size factor, and "DESeq2" takes the median ratio of a feature to its geometric
#'   mean across all samples
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
#' @param base the base used for the logarithm. Currently only supports
#'  "e" or "2" (can also specify the number 2).
#'
#' @return a list with two items:
#'  \itemize{
#'    \item{n_func}{the normalization function to be used with the
#'      'norm_fun_counts' and 'norm_fun_tpm' options in sleuth}
#'    \item{t_func}{the transformation function to be used with the
#'      'transform_fun_counts' and 'transform_fun_tpm' options in sleuth}
#'  }
#' 
#' @importFrom utils packageVersion compareVersion
#' @export
get_lr_functions <- function(type = "alr", denom_name = NULL,
                             lr_method = "both",
                             denom_method = "geomean",
                             impute_method = "multiplicative",
                             delta = NULL, impute_proportion = 0.65,
                             base = "e") {
  base <- as.character(base)
  base <- match.arg(base, c("e", "2"))
  denom_method <- match.arg(denom_method, c("geomean", "DESeq2"))
  impute_method <- match.arg(impute_method, c("multiplicative", "additive"))
  lr_method <- match.arg(lr_method, c("both", "transform"))
  type <- match.arg(type, c("alr", "ALR", "clr", "CLR", "iqlr", "IQLR"))
  type <- tolower(type)
  
  if (type == "alr" & is.null(denom_name)) {
    stop("you selected the 'ALR' transformation, but no ",
         "denominator was supplied. 'denom_name' is required ",
         "for the 'ALR' type.")
  }

  sleuth_version <- as.character(utils::packageVersion('sleuth'))
  sleuth_check <- utils::compareVersion(sleuth_version, '0.30.0')
  if (lr_method == 'both' && sleuth_check != 1) {
    msg <- paste0("The version of sleuth loaded is '", sleuth_version, "'. ",
                  "The lr_method 'both' only works with the API of sleuth ",
                  "in versions > 0.30.0. Setting the 'lr_method' to ",
                  "'transform'...")
    warning(msg)
    lr_method <- 'transform'
  }

  transform_func <- retrieve_transform_func(type = type, lr_method = lr_method,
                                            denom = denom_name,
                                            delta = delta,
                                            denom_method = denom_method,
                                            impute_proportion = impute_proportion,
                                            impute_method = impute_method,
                                            base = base)
  if (lr_method == "both") {
    message(">> ", Sys.time(), " - preparing sleuth object using the ",
            "sequential normalize and transform approach")
    norm_func <- retrieve_norm_func(type = type, denom = denom_name,
                                    delta = delta,
                                    denom_method = denom_method,
                                    impute_proportion = impute_proportion,
                                    impute_method = impute_method,
                                    base = base)
  } else {
    message(">> ", Sys.time(), " - preparing sleuth object using the ",
            "all-in-one transformation approach")
    norm_func <- norm_identity
  }

  return(list(n_func = norm_func, t_func = transform_func))
}

retrieve_transform_func <- function(type = "alr", lr_method = "both",
                                    denom = NULL, delta = 0.01,
                                    denom_method = "geomean",
                                    impute_proportion = 0.65,
                                    impute_method = "multiplicative",
                                    base = "e")
{
  e <- new.env()
  e$delta <- delta
  e$impute <- impute_proportion
  e$impute_method <- impute_method
  e$base <- base

  if (lr_method == "transform") {
    e$denom_method <- denom_method
    if (type == "alr") {
      e$denom <- denom
      e$fun <- function(matrix, sf = 1, denom_name = eval(e$denom)) {
        alr_transformation(matrix, denom_name = denom_name,
                           base = e$base, delta = e$delta,
                           denom_method = e$denom_method,
                           impute_method = e$impute_method,
                           impute_proportion = e$impute)
      }
    } else if (type == "iqlr") {
      e$denom <- "iqlr"
      e$fun <- function(matrix, sf = 1) {
        iqlr_transformation(matrix,
                            base = e$base, delta = e$delta,
                            denom_method = e$denom_method,
                            impute_method = e$impute_method,
                            impute_proportion = e$impute)
      }
    } else {
      e$denom <- "all"
      e$fun <- function(matrix, sf = 1) {
        clr_transformation(matrix,
                           base = e$base, delta = e$delta,
                           denom_method = e$denom_method,
                           impute_method = e$impute_method,
                           impute_proportion = e$impute)
      }
    }
  } else if (lr_method == "both") {
    e$fun <- function(mat, sf = 1, delta = e$delta, impute = e$impute,
                      impute_method = e$method, base = e$base) {
      logfunc <- switch(base, "e" = log, "2" = log2)
      mat <- suppressWarnings(
        impute_zeros(mat, delta = delta,
                     impute_proportion = impute_proportion,
                     method = impute_method)
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
  } else {
    stop("'lr_method' must be 'both' or 'transform'")
  }

  transform_func <- function(mat, sf) { e$fun(mat, sf = sf) }
  environment(transform_func) <- e

  return(transform_func)
}


retrieve_norm_func <- function(type = "alr", denom = NULL,
                               impute_method = "multiplicative",
                               denom_method = "geomean",
                               delta = NULL, impute_proportion = 0.65,
                               base = "e")
{
  n <- new.env()
  if (type != "alr") {
    denom <- NULL
  }

  n$denom_name <- denom
  n$type <- type
  n$method <- denom_method
  n$fun <- function(mat, denoms = n$denom_name, type = n$type, method = n$method) {
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

  norm_func
}
