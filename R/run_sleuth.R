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
#' @param aggregate_column, character indicating the column in
#'   \code{target_mapping} to be used for "gene-level" sleuth analysis.
#' @param num_cores, the number of cores to be used for sleuth analysis
#' @param lr_type, either "alr" or "clr" ("ALR" / "CLR" also accepted),
#'   indicating additive logratio or centered logratio transformation
#' @param denom_name, target ID names or index numbers of denominators;
#'   required for the ALR transformation.
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
                                  aggregate_column = NULL,
                                  num_cores = parallel::detectCores() - 2,
                                  lr_type = "alr", denom_name = NULL, ...)
{
  extra_opts <- list(...)
  if ("read_bootstrap_tpm" %in% names(extra_opts))
    read_bootstrap_tpm <- extra_opts$read_bootstrap_tpm
  else
    read_bootstrap_tpm <- TRUE

  if ("extra_bootstrap_summary" %in% names(extra_opts))
    extra_bootstrap_summary <- extra_opts$extra_bootstrap_summary
  else
    extra_bootstrap_summary <- TRUE

  if (is.null(extra_opts$impute_proportion))
    impute_proportion <- 0.65
  else impute_proportion <- extra_opts$impute_proportion

  # make the sleuth object using the PREP method,
  # which downloads the kallisto results and initializes the sleuth object
  # see ?sleuth::sleuth_prep for additional details
  transform_function <- get_lr_function(type = lr_type, denom_name = denom_name,
                                        delta = extra_opts$delta,
                                        impute_proportion = impute_proportion)
  sleuth.obj <- sleuth::sleuth_prep(sample_to_covariates, full_model,
                            target_mapping = target_mapping,
                            norm_fun_counts = norm_identity,
                            norm_fun_tpm = norm_identity,
                            aggregation_column = aggregate_column,
                            read_bootstrap_tpm = read_bootstrap_tpm,
                            extra_bootstrap_summary = extra_bootstrap_summary,
                            transformation_function = transform_function,
                            num_cores = num_cores, ...)
  # the default of sleuth_fit is to fit the 'full' model,
  # found in the 'full_model' variable above
  # see ?sleuth::sleuth_fit for more details
  sleuth.obj <- sleuth::sleuth_fit(sleuth.obj)
  # use Wald test to see significance of the chosen 'beta' on expression
  if (!is.null(beta))
    sleuth.obj <- sleuth::sleuth_wt(sleuth.obj, beta)
  message(">> ", Sys.time(), " - fitting null model and ",
          "performing likelihood ratio test")
  # use likelihood ratio test to look at the null_model versus the full_model
  sleuth.obj <- sleuth::sleuth_fit(sleuth.obj, formula = null_model,
                                   fit_name = "reduced")
  sleuth.obj <- sleuth::sleuth_lrt(sleuth.obj, "reduced", "full")
  sleuth.obj
}

