#' Make Compositional Sleuth Object
#'
#' This is a wrapper function for \code{sleuth_prep} that applies
#' the compositional data analysis approach. Many of the arguments
#' are ones that will be input into \code{sleuth_prep}.
#' 
#' @param sample_to_covariates, the sample_to_covariates matrix for sleuth
#' @param full_model, the full model for sleuth
#' @param target_mapping, the target mapping data frame for sleuth
#' @param beta, the beta you wish to use for the Wald test. If NULL,
#'   the Wald test will be skipped.
#' @param null_model, the null model to be the baseline for the LR test.
#'   default is the intercept only model (~1).
#' @param run_models boolean to see if the modeling step should be done.
#'   If \code{FALSE}, only sleuth_prep is done. Default is \code{TRUE}.
#' @param aggregate_column, character indicating the column in
#'   \code{target_mapping} to be used for "gene-level" sleuth analysis.
#' @param num_cores, the number of cores to be used for sleuth analysis
#' @param lr_type, either "alr" or "clr" ("ALR" / "CLR" also accepted),
#'   indicating additive logratio or centered logratio transformation
#' @param denom_name, target ID names or index numbers of denominators;
#'   required for the ALR transformation. If 'clr' or 'iqlr' is used, this overrides
#'   'lr_type' and the CLR / IQLR transformation will be used. If lr_type is 'clr' or 'iqlr',
#'   this argument will be ignored. If 'best' is used, then the internal function
#'   choose_best_denoms will be used to identify the feature(s) with the most consistent
#'   abundance across all samples in the experiment. See \link{choose_denom} for extra
#'   options.
#' @param which_var, must be "obs_tpm" or "obs_counts", to indicate
#'   whether sleuth should model TPMs or estimated counts, respectively
#' @param delta a number that is the imputed value during the transformation, to
#'  avoid zeros. If \code{NULL}, delta = impute_proportion * (minimum value in sample)
#' @param impute_proportion percentage of minimum value that
#'  becomes the imputed value. Only used if delta is \code{NULL} 
#' @param base the base used for the logarithm. Currently only supports
#'  "e" or "2" (can also specify the number 2).
#' @param ... extra options that will be passed on the sleuth.
#'   you can specify here whether \code{read_bootstrap_tpm} and
#'   \code{extra_bootstrap_summary} should be \code{FALSE} (default \code{TRUE}).
#' @return a sleuth object that has been prepped and fitted using the
#'   full and null models. It will also run the Wald test (if applicable)
#'   and the LR test
#'
#' @export
make_lr_sleuth_object <- function(sample_to_covariates, full_model = stats::formula('~condition'),
                                  target_mapping, beta, null_model = stats::formula('~1'),
                                  run_models = TRUE, aggregate_column = NULL,
                                  num_cores = parallel::detectCores() - 2,
                                  lr_type = "alr", denom_name = NULL, which_var = "obs_tpm",
                                  delta = NULL, impute_proportion = 0.65, base = "e", ...)
{
  stopifnot(which_var %in% c('obs_tpm', 'obs_counts'))
  best_denom_var <- ifelse(which_var == 'obs_tpm', 'tpm', 'est_counts')

  extra_opts <- list(...)
  if ("read_bootstrap_tpm" %in% names(extra_opts))
    read_bootstrap_tpm <- extra_opts$read_bootstrap_tpm
  else
    read_bootstrap_tpm <- TRUE

  if ("extra_bootstrap_summary" %in% names(extra_opts))
    extra_bootstrap_summary <- extra_opts$extra_bootstrap_summary
  else
    extra_bootstrap_summary <- TRUE

  if (!is.null(denom_name) && length(denom_name) == 1 && denom_name == 'best') {
    denom_name <- choose_denom(sample_info = sample_to_covariates,
                               target_mapping = target_mapping,
                               aggregation_column = aggregate_column,
                               num_cores = num_cores,
                               which_var = best_denom_var)
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

  # make the sleuth object using the PREP method,
  # which downloads the kallisto results and initializes the sleuth object
  # see ?sleuth::sleuth_prep for additional details
  transform_function <- get_lr_function(type = lr_type, denom_name = denom_name,
                                        delta = delta, impute_proportion = impute_proportion,
                                        base = base)
  sleuth.obj <- sleuth::sleuth_prep(sample_to_covariates, full_model,
                            target_mapping = target_mapping,
                            norm_fun_counts = norm_identity,
                            norm_fun_tpm = norm_identity,
                            aggregation_column = aggregate_column,
                            read_bootstrap_tpm = read_bootstrap_tpm,
                            extra_bootstrap_summary = extra_bootstrap_summary,
                            transform_fun_counts = transform_function,
                            transform_fun_tpm = transform_function,
                            num_cores = num_cores, ...)

  sleuth.obj <- clean_denom_names(sleuth.obj)

  # the default of sleuth_fit is to fit the 'full' model,
  # found in the 'full_model' variable above
  # see ?sleuth::sleuth_fit for more details
  if (run_models) {
    sleuth.obj <- sleuth::sleuth_fit(sleuth.obj, which_var = which_var)
    # use Wald test to see significance of the chosen 'beta' on expression
    if (!is.null(beta))
      sleuth.obj <- sleuth::sleuth_wt(sleuth.obj, beta)
    message(">> ", Sys.time(), " - fitting null model and ",
            "performing likelihood ratio test")
    # use likelihood ratio test to look at the null_model versus the full_model
    sleuth.obj <- sleuth::sleuth_fit(sleuth.obj, formula = null_model,
                                     fit_name = "reduced", which_var = which_var)
    sleuth.obj <- sleuth::sleuth_lrt(sleuth.obj, "reduced", "full")
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
