# sleuth-CN

Preprint available [here](https://www.biorxiv.org/content/10.1101/564955v1)! If you use this tool, any constructive feedback is greatly appreciated.

## Compositional Data Analysis approach to Sleuth

This is an extension of [`sleuth`](https://github.com/pachterlab/sleuth)
that treats the data as compositional data.
Applying the ideas of John Aitchison, as well as the ideas found in
the R packages `ALDEx2` and `compositions`, this performs a logratio
transformation of data created by `kallisto` or `salmon` to be 
analyzed by `sleuth`.

The "CN" stands for **Compositional Normalization**. Note that this approach is 
different from a formal compositional data analysis approach, which focuses its
attention purely on the relative abundances. This approach explicitly attempts
to recover the absolute abundances by anchoring the data on defined
"negative control" features (spike-ins, reference genes, etc). Even when such
features are not present in the data, it provides a more accurate interpretation
of the data (how all other features are changing *relative to* a defined set of features).

See the preprint for full details.

## Installation

To install, use devtools:

```
devtools::install_github('https://github.com/warrenmcg/sleuth-ALR')
```

## Usage

See the [vignette](https://github.com/warrenmcg/sleuth-ALR/blob/master/vignettes/sleuthALR.Rmd)
for full details!

Compared to running `sleuth_prep` and `sleuth_fit`, the only new arguments are the following:
+ `null_model`: specify both the full and null model, in order to run an LRT
+ `lr_type`: this specifies what kind of logratio transformation you want.
  use `alr` for additive logratio transformation, and use `clr` for centered
  logratio transformation (the geometric mean in each sample is the denominator).
+ `denom_name`: this must be specified if you choose `alr` for your `lr_type`.
  This is the target ID (or index number of the target ID) that you wish to be
  used as a 'reference gene'. If you specify more than one, the geometric mean
  between all of the selected reference genes will be used as the denominator.
+ additional arguments taken by `choose_denom` and `get_lr_functions`. See the help pages
for those functions for additional details.

## Options to handle zeros for ALR transformation

The ALR transformation involves the logarithm of a ratio, so a zero in either the
numerator or the denominator is a problem. The question of how to handle zeros remains
an area of active research within Compositional Data Analysis. For a review of this work,
please see [Chapter 4 of *Compositional Data Analysis*](http://onlinelibrary.wiley.com/doi/10.1002/9781119976462.ch4/summary).

sleuth-ALR takes the following approach:

1) **Essential Zeros**: If any feature is zero in all samples, or is otherwise filtered out by the sleuth
filtering step, it is considered an "essential zero". Biologically, the interpretation is
that this feature is likely in a region of heterochromatin and is silenced. These features
can be safely excluded from analysis.

2) **Rounded Zeros**: Any feature that has zero in at least one sample, but otherwise passes the sleuth
filtering step, it is considered a "rounded zero". Biologically, the interpretation is that this
feature is likely in a region where expression is occurring, but is below the limit of detection. The samples
with zeros are imputed using the "multiplicative strategy" described in Chapter 4 (link above).

Briefly, the features within a sample that are zero are replaced with `delta`. The non-zero features within that sample
are replaced by the following formula:

*old_value* = *x*   ==>   *new_value* = x \* (1 - (# of zero features) \* ***delta*** / *sum constraint*)

The `sum constraint` is the total for that composition (if TPMs, it would be 1 million; for counts, it would be the library size). 

Two additional arguments that can be specified are used in the ALR transformation:
+ `delta`: you can specify delta to be used for imputing (default: `NULL` and `impute_proportion` is used)
+ `impute_proportion`: if `delta` is `NULL`, the minimum non-zero value is taken to be the *detection limit*, and
`impute_proportion` is what proportion of this detection limit should be used as ***delta***. (default is 65%).

Chapter 4 from *Compositional Data Analysis* mentions the recommendation of 65%, hence the choice of default. However,
our own sensitivity analysis with knockout datasets indicates that a default value of ***delta*** = **0.5** minimizes
unstable variation in low-expression samples when there is a significant down-regulation of expression.
