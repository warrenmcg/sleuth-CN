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

#' Get Denom Name(s)
#' 
#' Function to retrieve the denominator(s) for a sleuth object
#' processed using compositional data analysis. If the analysis
#' was done using CLR (centered logratio), the denominator name is
#' all.
#'
#' @param obj the sleuth object
#' @return The name of the denominator(s) for the sleuth object
#' @export
get_denom_names <- function(obj) {
  stopifnot(is(obj, "sleuth"))
  name <- environment(obj$transform_fun)$denom_name
  if (is.null(name)) {
    stop("This sleuth object was not processed using compositional data analysis.",
         "\nThus its denominator does not exist.")
  }
  name
}

# function to clean the sleuth object to remove the denominator names
# from the processed data to prevent the denominator from affecting the modeling
# steps
clean_denom_names <- function(obj) {
  stopifnot(is(obj, 'sleuth'))
  denom_names <- get_denom_names(obj)
  if (is.null(denom_names) | denom_names == "all") {
    return(obj)
  } else {
    if(!is.null(obj$bs_summary$obs_counts)) {
      indices <- which(rownames(obj$bs_summary$obs_counts) %in% denom_names)
      obj$bs_summary$obs_counts <- obj$bs_summary$obs_counts[-indices, ]
      obj$bs_summary$sigma_q_sq <- obj$bs_summary$sigma_q_sq[-indices]
    }
    if(!is.null(obj$bs_summary$obs_tpm)) {
      indices <- which(rownames(obj$bs_summary$obs_tpm) %in% denom_names)
      obj$bs_summary$obs_tpm <- obj$bs_summary$obs_tpm[-indices, ]
      obj$bs_summary$sigma_q_sq_tpm <- obj$bs_summary$sigma_q_sq_tpm[-indices]
    }
    if(length(obj$bs_quants) > 0 & length(obj$bs_quants[[1]]) > 0) {
      obj$bs_quants <- lapply(obj$bs_quants, function(sample) {
        lapply(sample, function(summary) {
          indices <- which(rownames(summary) %in% denom_names)
          summary[-indices, ]
        })
      })
    }
    return(obj)
  }
}

#' Choose denominator
#'
#' This function selects the best denominator for an ALR
#' transformation based on which feature has the most consistent
#' proportion across all samples. Currently, only "coefficient of variation"
#' (method "cov") is implemented. Users can supply a "sleuth" object, or
#' alternatively can supply a sample_to_covariates table. It runs a truncated
#' form of sleuth_prep to get raw counts and TPMs to calculate the metric.
#'
#' @param obj a sleuth object. This or sample_info must be specified. If both
#'   are specified, the sleuth object will take precedence.
#' @param sample_info a sample_to_covariates table that must match sleuth API.
#'   see ?sleuth::sleuth_prep for more information
#' @param target_mapping a table that matches each target_id to gene identifiers
#'   or other metadata. If aggregation_column is specified, this table must be
#'   provided, and it must have aggregation_column as one of its named columns.
#' @param aggregation_column If you plan to examine gene-level information,
#'   you must provide an aggregation_column to aggregate transcripts to the gene level.
#'   You must provide a target_mapping table with a column whose name is this argument.
#' @param num_cores the number of cores that parallel should use to process samples
#' @param which_var must be one of "est_counts", "tpm", or "scaled_reads_per_base" (for gene-level counts).
#' @param method the metric used to select the denominator. Currently only "cov" is implemented.
#'
#' @return a character vector with one or more denominators for use. If multiple
#'   features have ties for the metric of choice, currently all of them are returned so
#'   that their geometric mean can be used as the denominator.
#' @export
choose_denom <- function(obj = NULL, sample_info = NULL, target_mapping = NULL,
                         aggregation_column = NULL, num_cores = 1,
                         which_var = "est_counts", method = "cov") {
  if (is.null(obj)) {
    if (is.null(sample_info)) {
      stop("You must provide either a sleuth object to 'obj', ",
           "or a sample_to_covariates data frame to 'sample_info'.")
    } else if (!is.null(aggregation_column) & is.null(target_mapping)) {
      stop("If an aggregation_column is provided, you must also ",
           "provide target_mapping to do the aggregation on.")
    }
    message("Preparing the sleuth object to select the best denominator")
    obj <- suppressMessages(sleuth::sleuth_prep(sample_info,
                                                target_mapping = target_mapping,
                                                aggregation_column = aggregation_column,
                                                norm_fun_counts = norm_identity,
                                                norm_fun_tpm = norm_identity,
                                                max_bootstrap = 1, num_cores = 1,
                                                full_model = formula("~1")))
  } else {
    if(!is(obj, "sleuth")) {
      stop("You must provide either a sleuth object to 'obj', ",
           "or a sample_to_covariates data frame to 'sample_info'.")
    } else if (!is.null(sample_info)) {
      warning("You provided both a sleuth object and a ",
              "sample_to_covariates data frame. The s2c data frame will be ignored.")
    }
  }

  stopifnot(which_var %in% c("est_counts", "scaled_reads_per_base", "tpm"))
  which_var <- sleuth:::check_quant_mode(obj, which_var)

  if (method == "cov") {
    message("Calculating the coefficient of variation of all targets")
    if (is.null(aggregation_column)) {
      mat <- sleuth:::spread_abundance_by(obj$obs_raw, which_var)
    } else {
      mat <- sleuth:::spread_abundance_by(obj$obs_norm, which_var)
    }
    cov <- apply(mat, 1, function(x) sd(x, na.rm = T) / mean(x, na.rm = T))
    min_cov <- min(cov, na.rm = T)
    min_index <- which(cov == min_cov)
    denom_name <- names(cov)[min_index]
  } else {
    stop("Currently, 'cov' is the only implemented method. ",
         "Contact the developer if you are interested in another method.")
  }
  denom_name
}
