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
#'   you must provide an aggregation_column to aggregate transcripts to the gene
#'   level. You must provide a target_mapping table with a column whose name is
#'   this argument or a sleuth object with that table.
#' @param num_cores the number of cores that parallel should use to process samples
#' @param num_denoms the number of features to select for normalization. The default is one.
#' @param which_var must be one of "est_counts", "tpm", or "scaled_reads_per_base" (for gene-level counts).
#' @param min_value the minimum threshold for the mean of 'which_var' for the candidate denominator.
#' @param method the metric used to select the denominator. Currently only "cov" is implemented.
#' @param filter_length boolean to filter possible denominators by length. This requires that 'length'
#'   is a column in the target_mapping. If \code{TRUE}, then all features with a length less
#'   than 300 bases is excluded from consideration. This is recommended when modeling TPMs,
#'   and is thus \code{TRUE} by default since TPMs are modeled by default.
#'
#' @return a character vector with one or more denominators for use. If multiple
#'   features have ties for the metric of choice, currently all of them are returned so
#'   that their geometric mean can be used as the denominator.
#' @export
choose_denom <- function(obj = NULL, sample_info = NULL, target_mapping = NULL,
                         aggregation_column = NULL, num_cores = 1, num_denoms = 1,
                         which_var = "tpm", min_value = 5, method = "cov",
                         filter_length = TRUE) {
  stopifnot(which_var %in% c("est_counts", "scaled_reads_per_base", "tpm"))
  if (is.null(obj)) {
    if (is.null(sample_info)) {
      stop("You must provide either a sleuth object to 'obj', ",
           "or a sample_to_covariates data frame to 'sample_info'.")
    } else if (!is.null(aggregation_column) & is.null(target_mapping)) {
      stop("If an aggregation_column is provided, you must also ",
           "provide target_mapping to do the aggregation on.")
    }
    message("Preparing the sleuth object to select the best denominator")
    if (!is.null(aggregation_column)) {
      obj <- suppressMessages(sleuth::sleuth_prep(sample_info,
                                                  target_mapping = target_mapping,
                                                  aggregation_column = aggregation_column,
                                                  norm_fun_counts = norm_identity,
                                                  norm_fun_tpm = norm_identity,
                                                  max_bootstrap = 2, num_cores = 1,
                                                  full_model = formula("~1")))
    } else {
      obj <- suppressMessages(sleuth::sleuth_prep(sample_info,
                                                  target_mapping = target_mapping,
                                                  num_cores = 1,
                                                  normalize = FALSE))
    }
  } else {
    if(!is(obj, "sleuth")) {
      stop("You must provide either a sleuth object to 'obj', ",
           "or a sample_to_covariates data frame to 'sample_info'.")
    } else if (!is.null(sample_info)) {
      warning("You provided both a sleuth object and a ",
              "sample_to_covariates data frame. The s2c data frame will be ignored.")
    }
  }

  if (!is.null(aggregation_column) && which_var == "est_counts") {
    which_var <- "scaled_reads_per_base"
  } else if (is.null(aggregation_column) && which_var == "scaled_reads_per_base") {
    which_var <- "est_counts"
  }

  if (filter_length) {
    message("filtering by length")
    if(!("length" %in% names(target_mapping))) {
      stop("The 'target_mapping' table is missing a 'length' column. This is required to filter features by length.")
    }
    if(!is.null(aggregation_column)) {
      transcript_lengths <- data.table::as.data.table(dplyr::select_(target_mapping, aggregation_column, "length"))
      length_bool <- transcript_lengths[, list(bool = all(length >= 300)), by = list(eval(parse(text = aggregation_column)))]
      length_bool_ids <- length_bool[[aggregation_column]][length_bool$bool]
    } else {
      transcript_lengths <- dplyr::select(target_mapping, target_id, length)
      length_bool <- transcript_lengths$length >= 300
      length_bool_ids <- transcript_lengths$target_id[length_bool]
    }
  }

  if (method == "cov") {
    message("Calculating the coefficient of variation of all targets")
    if (is.null(aggregation_column)) {
      mat <- sleuth::sleuth_to_matrix(obj, 'obs_raw', which_var)
    } else {
      mat <- sleuth::sleuth_to_matrix(obj, 'obs_norm', which_var)
    }
    # only have a matrix containing filtered values, to remove low expressing transcripts
    # which often have low coefficients of variation
    # gene-level, using the normalized matrix, is already filtered
    allNonZero <- !matrixStats::rowAnys(mat, value = 0)
    mat <- mat[obj$filter_bool & allNonZero, ]

    if (filter_length) mat <- mat[row.names(mat) %in% length_bool_ids, ]
    mean_vals <- rowMeans(mat, na.rm = T)
    cov <- matrixStats::rowSds(mat, na.rm = T) / rowMeans(mat, na.rm = T)

    if (num_denoms == 1) {
      min_cov <- min(cov[which(mean_vals >= min_value & cov > 0)], na.rm = T)
      min_index <- which(cov == min_cov)
    } else {
      filt_cov <- cov[which(mean_vals >= min_value & cov > 0)]
      ordered_cov <- filt_cov[order(filt_cov)]
      min_cov <- ordered_cov[1:num_denoms]
      min_index <- which(cov %in% min_cov)
    }
    denom_name <- names(cov)[min_index]
  } else {
    stop("Currently, 'cov' is the only implemented method. ",
         "Contact the developer if you are interested in another method.")
  }
  denom_name
}
