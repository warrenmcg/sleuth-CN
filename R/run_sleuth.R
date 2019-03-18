#' Make Compositional Sleuth Object
#'
#' This is a wrapper function for the sleuth pipeline that applies
#' the compositional normalization approach. Many of the arguments
#' are ones that will be input into different parts of the sleuth
#' pipeline: \code{sleuth_prep}, \code{sleuth_fit}, \code{sleuth_wt},
#' and \code{sleuth_lrt}.
#' 
#' @param sample_to_covariates, the sample_to_covariates \code{data.frame} for
#'   sleuth
#' @param denom_name, target ID names or index numbers of features to be used
#'   for the denominator when using compositional analysis; this argument is
#'   required for the ALR transformation. Using 'clr' or 'iqlr' for 'lr_type' 
#'   overrides this argument (and the CLR / IQLR transformation will be used).
#'   If 'best' is used, then the internal function \code{choose_best_denoms}
#'   will be used to identify the feature(s) with the most consistent abundance
#'   across all samples in the experiment. See \link{choose_denom} for extra
#'   options.
#' @param lr_type, either "alr", "clr", or "iqlr" ("ALR" / "CLR" / "IQLR" are
#'   also accepted), indicating additive, centered, or interquartile logratio
#'   transformation
#' @param full_model, the full model for sleuth
#' @param beta, the beta you wish to use for the Wald test. If \code{NULL}
#'   (the default), the Wald test will be skipped.
#' @param null_model, the null model to be the baseline for the LR test.
#'   If \code{NULL} (the default), this step is skipped.
#' @param run_models boolean to see if the modeling step should be done.
#'   If \code{FALSE}, the default, only sleuth_prep is done.
#' @param ... extra options that will tweak the analysis, specifically for
#'  \link{\code{get_lr_functions}}, \code{sleuth_prep}, and \code{sleuth_fit}.
#'  for details on which options can be specified for \code{sleuth_prep} and
#'  \code{sleuth_fit}, please see ?get_lr_functions, ?sleuth::sleuth_prep or
#'  ?sleuth::sleuth_fit for details. Note that for \code{sleuth_prep},
#'  \code{read_bootstrap_tpm} and \code{extra_bootstrap_summary} are \code{TRUE}
#'  by default to allow for modeling estimated or TPMs downstream.
#' @return a sleuth object that has been prepped using the compositional
#'   analysis for the normalization and transformation steps, and fitted
#'   using the full model (if \code{run_models} is \code{TRUE}) and null model
#'   (if specified). It will also run the Wald test (if \code{beta} is
#'   specified) and the LR test (if applicable).
#'
#' @export
make_lr_sleuth_object <- function(sample_to_covariates, denom_name,
                                  lr_type = "alr", full_model = NULL,
                                  beta = NULL, null_model = NULL,
                                  run_models = FALSE, ...)
{
  opts_list <- process_extra_opts(...)
  prep_opts <- opts_list$prep_opts
  fit_opts <- opts_list$fit_opts
  alr_opts <- opts_list$alr_opts
  choose_opts <- opts_list$choose_opts

  if (run_models && !is.null(extra_opts$beta)) {
    test_design <- model.matrix(full_model, data = sample_to_covariates)
    if (!beta %in% colnames(test_design)) {
      stop("It does not seem that the specified 'beta' is found among the ",
           "possible betas for your specified 'sample_to_covariates' and ",
           "'full_model'. Here are the possible betas: ",
           paste(colnames(test_design), collapse = ", "))
    }
  } else if (run_models && !is.null(full_model)) {
    stop("If 'run_models' is TRUE, then 'full_model' must be specified")
  }

  if (!is.null(denom_name) && length(denom_name) == 1 && denom_name == 'best') {
    opts <- list(sample_info = sample_to_covariates,
                 target_mapping = prep_opts$target_mapping,
                 aggregation_column = prep_opts$aggregation_column,
                 num_cores = prep_opts$num_cores)
    choose_opts <- c(opts, choose_opts)
    denom_name <- do.call(choose_denom, choose_opts)
  } else if (!is.null(denom_name) && tolower(denom_name) == 'clr') {
    if(tolower(lr_type) != 'clr') {
      message("'denom_name' is 'clr', but 'lr_type' is not 'clr'. The lr_type will be overriden and is now 'clr'")
      lr_type <- 'clr'
    }
  } else if (!is.null(denom_name) && tolower(denom_name) == 'iqlr') {
    if(tolower(lr_type) != 'iqlr') {
      message("'denom_name' is 'iqlr', but 'lr_type' is not 'iqlr'. The lr_type will be overriden and is now 'iqlr'")
      lr_type <- 'iqlr'
    }
  }

  # create the ALR normalization and transformation functions
  # see `?get_lr_functions` for additional details
  opts <- list(type = lr_type, denom_name = denom_name)
  alr_opts <- c(opts, alr_opts)

  func_list <- do.call(get_lr_functions, args = alr_opts)
  norm_func <- func_list$n_func
  transform_func <- func_list$t_func

  # make the sleuth object using the `sleuth_prep` function,
  # which downloads the kallisto results and initializes the sleuth object
  # see `?sleuth::sleuth_prep` for additional details
  sleuth_opts <- list(sample_to_covariates = sample_to_covariates,
                      full_model = full_model,
                      norm_fun_counts = norm_func,
                      norm_fun_tpm = norm_func,
                      transform_fun_counts = transform_func,
                      transform_fun_tpm = transform_func)

  prep_opts <- c(sleuth_opts, prep_opts)
  sleuth.obj <- do.call(sleuth::sleuth_prep, args = prep_opts)

  # this line of code makes sure that, if a single feature is used as a
  # denominator, it is disqualified from the fitting steps, since it
  # is no longer dependent.
  sleuth.obj <- clean_denom_names(sleuth.obj, lr_method = alr_opts$lr_method)

  # the default of sleuth_fit is to fit the 'full' model,
  # found in the 'full_model' variable above
  # see ?sleuth::sleuth_fit for more details
  if (run_models) {
    message(">> ", Sys.time(), " - fitting full model")
    full_opts <- list(obj = sleuth.obj)
    fit_opts <- c(full_opts, fit_opts)
    sleuth.obj <- do.call(sleuth::sleuth_fit, args = fit_opts)
    # use Wald test to see significance of the chosen 'beta' on expression
    if (!is.null(beta)) {
      sleuth.obj <- sleuth::sleuth_wt(sleuth.obj, beta)
    } else {
      message("no 'beta' specified; skipping Wald test for the full model")
    }
    # use likelihood ratio test to look at the null_model versus the full_model
    if (!is.null(null_model)) {
      message(">> ", Sys.time(), " - fitting null model and ",
              "performing likelihood ratio test")
      null_opts <- list(formula = null_model, fit_name = "reduced")
      fit_opts <- c(fit_opts, null_opts)
      sleuth.obj <- do.call(sleuth::sleuth_fit, fit_opts)
      sleuth.obj <- sleuth::sleuth_lrt(sleuth.obj, "reduced", "full")
    } else {
      message("no 'null_model' specified; skipping null model fitting and ",
              "likelihood ratio test")
    }
  } else {
    message("'run_models' is FALSE; skipping fitting and testing steps")
  }
  sleuth.obj
}

#' Return sleuth-ALR Results Table
#'
#' This is a wrapper function for \code{sleuth_results} that properly
#' handles the unique features of sleuth-ALR-transformed data.
#'
#' @param obj a \code{sleuth} object
#' @param test a character string denoting the test to extract. Possible tests can be found by using \code{models(obj)}.
#' @param test_type 'wt' for Wald test or 'lrt' for Likelihood Ratio test.
#' @param which_model a character string denoting the model. If extracting a wald test, use the model name.
#'   Not used if extracting a likelihood ratio test.
#' @param show_all if \code{TRUE} will show all transcripts (not only the ones
#' passing filters). The transcripts that do not pass filters will have
#' \code{NA} values in most columns.
#' @param pval_aggregate if \code{TRUE} and both \code{target_mapping} and \code{aggregation_column} were provided,
#' to \code{sleuth_prep}, use lancaster's method to aggregate p-values by the \code{aggregation_column}.
#' @param weight_func if \code{pval_aggregate} is \code{TRUE}, then this is used to weight the p-values for
#'   lancaster's method. This must be either the string 'best' or it must be a function that takes the
#'   observed means of the transcripts as the only defined argument.
#'
#' @details For the transcript-level analysis, this produces the same results as the default \code{sleuth_results}.
#'   However, using the default function for p-value aggregation is incompatible with the standard sleuth-ALR
#'   transformation. Sleuth-ALR logratios typically include negative values (any feature that is less abundant than the
#'   chosen 'reference feature(s)' will yield negative logratios), and negative values are not allowed for the
#'   lancaster method.
#'
#'   This method works around this problem by specifying a weighting function that is compatible with the logratios
#'   and with the lancaster method. The default is to specify the string 'best', which uses an internal function to
#'   determine how to exponentiate the logratios to get the ratios, using whatever base was used for the transformation.
#'
#' @return a \code{data.frame} with the same specification as found in \code{sleuth_results}. See \code{sleuth_results} for
#'   details.
#' @export
sleuth_alr_results <- function(obj, test, test_type = "wt", which_model = "full",
                               show_all = TRUE, pval_aggregate = obj$pval_aggregate,
                               weight_func = 'best') {
  if (pval_aggregate) {
    if (weight_func == 'best') {
      weight_func <- get_alr_weight(obj)
    } else if (!is(weight_func, 'function')) {
      stop("'weight_func' must either be a function or the string 'best'.")
    }

    if (test_type == 'lrt') {
      res <- get_alr_test(obj, test, type = 'lrt')
    } else {
      res <- get_alr_test(obj, test, type = 'wt', model = which_model)
    }

    test_values <- weight_func(res$mean_obs)
    if (any(test_values < 0)) {
      stop("The selected 'weight_func' yielded negative values. This should not happen, so check ",
           "to make sure the correct function was used.")
    } else if (any(is.na(test_values))) {
      warning("The selected 'weight_func' yielded NA / NaN values. This should not typically happen, ",
              "so check to make sure the correct function was used. Continuing with retrieving the results anyway.")
    }
  }

  sleuth_results(obj, test = test, test_type = test_type, which_model = which_model,
                 show_all = show_all, pval_aggregate = pval_aggregate,
                 weight_func = weight_func)
}
