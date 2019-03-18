#' Process extra options for main sleuth-ALR function
#'
#' This internal function processes all of the extra arguments for various
#' steps of the sleuth-ALR pipeline.
#' 
#' @param ... extra arguments passed from \code{\link{make_lr_sleuth_object}}
#'
#' @details Four sets of arguments will be processed: sleuth_prep (prep_opts),
#'   sleuth_fit (fit_opts), get_lr_functions (alr_opts), and choose_denom
#'   (choose_opts). See the documentation for those functions if you want a
#'   list of arguments to use and how they are used.
process_extra_opts <- function(...) {

  extra_opts <- list(...)

  prep_opts <- list()
  fit_opts <- list()
  alr_opts <- list()
  choose_opts <- list()

  sleuth_version <- as.character(utils::packageVersion('sleuth'))
  sleuth_check <- utils::compareVersion(sleuth_version, '0.30.0')

  ##########################
  ### sleuth fit options ###
  ##########################

  # The following are extra arguments that can be specified in `sleuth_fit`
  # If the user does not supply them, then the defaults used by `sleuth_fit`
  # are specified
  if ("which_var" %in% names(extra_opts)) {
    fit_opts$which_var <- match.arg(extra_opts$which_var,
                                    c("obs_counts", "obs_tpm"))
    extra_opts$which_var <- NULL
  } else {
    fit_opts$which_var <- "obs_tpm"
  }

  if ("shrink_fun" %in% names(extra_opts)) {
    if (sleuth_check != 1) {
      msg <- paste0("The version of sleuth loaded is '", sleuth_version, "'. ",
                    "The 'shrink_fun' argument only works with the API of ",
                    "sleuth in versions > 0.30.0. Ignoring this option...")
      warning(msg)
    } else {
      if (!is.function(extra_opts$shrink_fun)) {
        stop("'shrink_fun' must be a function. Please see ?sleuth::sleuth_fit ",
             "for more details.")
      }
      fit_opts$shrink_fun <- extra_opts$shrink_fun
    }
  } else if (sleuth_check == 1) {
    fit_opts$shrink_fun <- sleuth::basic_shrink_fun
  }

  if ("n_bins" %in% names(extra_opts)) {
    fit_opts$n_bins <- extra_opts$n_bins
  } else {
    fit_opts$n_bins <- 100
  }

  if ("lwr" %in% names(extra_opts)) {
    lwr <- extra_opts$lwr
    if(!is.numeric(lwr) || lwr <= 0 || lwr >=1) {
      stop("'lwr' must between 0 and 1 (exclusive).")
    }
    fit_opts$lwr <- lwr
  } else {
    fit_opts$lwr <- 0.25
  }

  if ("upr" %in% names(extra_opts)) {
    upr <- extra_opts$upr
    if(!is.numeric(upr) || upr <= 0 || upr >=1) {
      stop("'upr' must between 0 and 1 (exclusive).")
    }
    fit_opts$upr <- upr
  } else {
    fit_opts$upr <- 0.75
  }

  if (fit_opts$upr <= fit_opts$lwr) {
    stop("'upr' must be larger than 'lwr'")
  }
  #################################
  ### end of sleuth fit options ###
  #################################

  ###################
  ### ALR options ###
  ###################
  if ("lr_method" %in% names(extra_opts)) {
    alr_opts$lr_method <- match.arg(extra_opts$lr_method,
                                    c("transform", "both"))
    extra_opts$lr_method <- NULL
  } else {
    alr_opts$lr_method <- "both"
  }

  if ("denom_method" %in% names(extra_opts)) {
    alr_opts$denom_method <- match.arg(extra_opts$denom_method,
                                       c("geomean", "DESeq2"))
    extra_opts$denom_method <- NULL
  } else {
    alr_opts$denom_method <- "geomean"
  }

  if ("impute_method" %in% names(extra_opts)) {
    alr_opts$impute_method <- match.arg(extra_opts$impute_method,
                                        c("multiplicative", "additive"))
    extra_opts$impute_method <- NULL
  } else {
    alr_opts$impute_method <- "multiplicative"
  }

  if ("delta" %in% names(extra_opts)) {
    stopifnot(is.numeric(extra_opts$delta))
    alr_opts$delta <- extra_opts$delta
    extra_opts$delta <- NULL
  } else {
    alr_opts$delta <- ifelse(fit_opts$which_var == "obs_tpm", 0.01, 0.5)
  }

  if ("impute_proportion" %in% names(extra_opts)) {
    if (!is.null(extra_opts$impute_proportion)) {
      stopifnot(is.numeric(extra_opts$impute_proportion))
      stopifnot(extra_opts$impute_proportion < 1 &&
                extra_opts$impute_proportion > 0)
    }
    alr_opts$impute_proportion <- extra_opts$impute_proportion
    extra_opts$impute_proportion <- NULL
  } else {
    alr_opts$impute_proportion <- 0.65
  }

  if ("base" %in% names(extra_opts)) {
    base <- as.character(extra_opts$base)
    alr_opts$base <- match.arg(base, c("e", "2"))
  } else {
    alr_opts$base <- "e"
  }

  ##########################
  ### end of ALR options ###
  ##########################

  ############################
  ### choose denom options ###
  ############################

  denom_var <- ifelse(fit_opts$which_var == 'obs_tpm', 'tpm', 'est_counts')
  if (!is.null(extra_opts$gene_mode) && denom_var == "est_counts") {
    denom_var <- "scaled_reads_per_base"
  }
  choose_opts$denom_var <- denom_var

  if ("num_denoms" %in% names(extra_opts)) {
    if (!is.numeric(extra_opts$num_denoms) || extra_opts$num_denoms <= 0) {
      stop("'num_denoms' must be a number that is 1 or more.")
    }
    choose_opts$num_denoms <- as.integer(extra_opts$num_denoms)
    extra_opts$num_denoms <- NULL
  } else {
    choose_opts$num_denoms <- 1
  }

  if ("min_value" %in% names(extra_opts)) {
    if (!is.numeric(extra_opts$min_value) || extra_opts$min_value <= 0) {
      stop("'min_value' must be a number that is 1 or more.")
    }
    choose_opts$min_value <- extra_opts$min_value 
    extra_opts$min_value <- NULL
  } else {
    choose_opts$min_value <- 5
  }

  if ("choose_method" %in% names(extra_opts)) {
    if (!is.character(extra_opts$choose_method) ||
          length(extra_opts$choose_method) != 1) {
      stop("'choose_method' must be a character vector with one value.")
    }
    choose_opts$choose_method <- extra_opts$choose_method
    extra_opts$choose_method <- NULL
  } else {
    choose_opts$choose_method <- "cov"
  }

  if ("filter_length" %in% names(extra_opts)) {
    if (!is.logical(extra_opts$filter_length) ||
          length(extra_opts$filter_length) != 1) {
      stop("'filter_length' must be a number that is 1 or mor.")
    }
    choose_opts$filter_length <- match.arg(extra_opts$filter_length,
                                           c(TRUE, FALSE))
    extra_opts$filter_length <- NULL
  } else {
    choose_opts$filter_length <- TRUE
  }

  ###################################
  ### end of choose denom options ###
  ###################################

  ###########################
  ### sleuth prep options ###
  ###########################
  if ("read_bootstrap_tpm" %in% names(extra_opts)) {
    prep_opts$read_bootstrap_tpm <- match.arg(extra_opts$read_bootstrap_tpm,
                                              c(TRUE, FALSE))
    extra_opts$read_bootstrap_tpm <- NULL
  } else {
    prep_opts$read_bootstrap_tpm <- TRUE
  }

  if ("extra_bootstrap_summary" %in% names(extra_opts)) {
    prep_opts$extra_bootstrap_summary <- match.arg(extra_opts$extra_bootstrap_summary,
                                                   c(TRUE, FALSE))
    extra_opts$extra_bootstrap_summary <- NULL
  } else {
    prep_opts$extra_bootstrap_summary <- TRUE
  }

  if ("num_cores" %in% names(extra_opts)) {
    prep_opts$num_cores <- extra_opts$num_cores
    extra_opts$num_cores <- NULL
  } else {
    prep_opts$num_cores <- 1
  }

  if ("gene_mode" %in% names(extra_opts)) {
    if (!is.logical(extra_opts$gene_mode) ||
          is.na(extra_opts$gene_mode)) {
      stop("'gene_mode' must be 'TRUE' or 'FALSE'")
    }
    prep_opts$gene_mode <- extra_opts$gene_mode
    extra_opts$gene_mode <- NULL
  } else {
    prep_opts$gene_mode <- FALSE
  }

  if ("target_mapping" %in% names(extra_opts)) {
    stopifnot(methods::is(extra_opts$target_mapping, "data.frame"))
    prep_opts$target_mapping <- extra_opts$target_mapping
    extra_opts$target_mapping <- NULL
  } else {
    prep_opts$target_mapping <- NULL
  }

  if ("aggregation_column" %in% names(extra_opts)) {
    agg_col <- extra_opts$aggregation_column
    if (is.null(prep_opts$target_mapping)) {
      stop("'target_mapping' is required if 'aggregation_column' is specified.")
    } else if (!agg_col %in% colnames(prep_opts$target_mapping)) {
      stop("'", agg_col, "' was not found in the column names for ",
           "'target_mapping'")
    }
    prep_opts$aggregation_column <- extra_opts$aggregation_column
    extra_opts$aggregation_column <- NULL
  } else {
    prep_opts$aggregation_column <- NULL
  }

  if (any(grepl("^norm_fun", names(extra_opts)))) {
    warning("a normalization function was specified as an extra option, ",
            "but this is produced interally using 'get_lr_functions'. ",
            "Ignoring ...")
    bad_opts <- grepl("norm_fun", names(extra_opts))
    extra_opts[bad_opts] <- NULL
  }

  if (any(grepl("transform_fun", names(extra_opts)))) {
    warning("a transformation function was specified as an extra option, ",
            "but this is produced internally using 'get_lr_functions'. ",
            "Ignoring ...")
    bad_opts <- grepl("transform_fun", names(extra_opts))
    extra_opts[bad_opts] <- NULL
  }

  ## any other options are intended to be for sleuth prep
  prep_opts <- c(prep_opts, extra_opts)
  ##################################
  ### end of sleuth prep options ###
  ##################################

  list(prep_opts = prep_opts, fit_opts = fit_opts, alr_opts = alr_opts, choose_opts = choose_opts)
}
