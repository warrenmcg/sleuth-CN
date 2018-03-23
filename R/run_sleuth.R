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
#' @param ... extra options that will be passed on the sleuth.
#'   you can specify here whether \code{read_bootstrap_tpm} and
#'   \code{extra_bootstrap_summary} should be \code{FALSE} (default \code{TRUE}).
#'   you can also specify here 'delta' and 'impute_proportion' for
#'   the transformation functions.
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
                                  delta = NULL, impute_proportion = 0.65, ...)
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

  if (!is.null(denom_name) && denom_name == 'best') {
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
                                        delta = delta, impute_proportion = impute_proportion)
  sleuth.obj <- sleuth::sleuth_prep(sample_to_covariates, full_model,
                            target_mapping = target_mapping,
                            norm_fun_counts = norm_identity,
                            norm_fun_tpm = norm_identity,
                            aggregation_column = aggregate_column,
                            read_bootstrap_tpm = read_bootstrap_tpm,
                            extra_bootstrap_summary = extra_bootstrap_summary,
                            transformation_function = transform_function,
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
