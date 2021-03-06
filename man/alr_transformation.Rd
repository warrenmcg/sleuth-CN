% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alr.R
\name{alr_transformation}
\alias{alr_transformation}
\title{Calculate ALR}
\usage{
alr_transformation(mat, denom_name, base = "e", remove_zeros = FALSE,
  denom_method = "geomean", impute_method = "multiplicative",
  delta = NULL, impute_proportion = 0.65)
}
\arguments{
\item{mat}{an D x M matrix of D target IDs and
M samples}

\item{denom_name}{a character vector of target IDs or
an integer vector of row numbers. If there is
more than one, the geometric mean of the
denominator values within one sample is used
as the denominator.}

\item{base}{what should the base of the logarithm be?
currently only supports base "e" and base 2.}

\item{remove_zeros}{boolean to see if this function
should remove essential zeros (features with zeros in
all samples). The default is \code{FALSE} to be
compatible with sleuth, as its default filter removes
essential zeros.}

\item{denom_method}{either 'geomean' or 'DESeq2' to
use either the geometric mean of the features as the
denominator, or the DESeq2-style size factors
as the denominator. If there is only one feature,
'geomean' uses the feature's value itself, whereas
'DESeq2' will use the ratio of the feature to 
the geometric mean across samples.}

\item{impute_method}{which method to use for imputing zeros.
'multiplicative' (default) sets all values smaller than
a imputation value 'delta' (determined by delta or
impute_proportion) to that imputation value, and reduces
all other values by the amount X * (1 - delta*num_zero_values /
sum_constraint). 'additive' is similar to most other tools, and
just adds the imputation value to all entries ('delta' must
be specified)}

\item{delta}{a number that is the imputed value. If \code{NULL},
delta = impute_proportion * (minimum value in sample)}

\item{impute_proportion}{percentage of minimum value that
becomes the imputed value. Only used if delta is \code{NULL}}
}
\value{
(D - n - z) x M matrix of ALR-transformed
  values, with n equal to the number of denominator values
  and z are the number of rows with essential zeros.
}
\description{
This function calculates the additive logratio
transformation. This function uses calculate_alr
internally after processing the data to handle
zeros. Essential zeros (i.e. target IDs with 
zeros in all samples) can be optionally excluded, and
all remaining zeros, including rounded zeros
(i.e. target IDs with at least one zero value) are
always imputed using the hidden impute_zeros function.
}
