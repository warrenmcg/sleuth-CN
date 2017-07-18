# sleuth.comp

## Compositional Data Analysis approach to Sleuth

This is an extension of [`sleuth`](https://github.com/pachterlab/sleuth)
that treats the data as compositional data.
Applying the ideas of John Aitchison, as well as the ideas found in
the R packages `ALDEx2` and `compositions`, this performs a logratio
transformation of data created by `kallisto` or `salmon` to be 
analyzed by `sleuth`.

## Installation

To install, use devtools:

```
devtools::install_github('https://github.com/pimentel/sleuth-lr')
```

## Usage

Using this package is very easy. Just use the wrapper function
`make_lr_sleuth_object`:

```
library(sleuth.comp)
make_lr_sleuth_object(sample_to_covariates, full_model = stats::formula('~condition'),
                      target_mapping, beta, null_model = stats::formula('~1'),
                      aggregate_column = NULL,
                      num_cores = parallel::detectCores() - 2,
                      lr_type = "alr", denom_name = NULL, ...)
```

Compared to running `sleuth_prep`, the only new arguments are the following:
+ `null_model`: specify both the full and null model, in order to run an LRT
+ `lr_type`: this specifies what kind of logratio transformation you want.
  use `alr` for additive logratio transformation, and use `clr` for centered
  logratio transformation (the geometric mean in each sample is the denominator).
+ `denom_name`: this must be specified if you choose `alr` for your `lr_type`.
  This is the target ID (or index number of the target ID) that you wish to be
  used as a 'reference gene'. If you specify more than one, the geometric mean
  between all of the selected reference genes will be used as the denominator.

### To-dos:

+ [ ] Finish the script to run gene-level "ground truth" comparison
+ [ ] Run a sensitivity analysis to choice of imputation value for "rounded zeros"
(default right now is 0.5, which is really high compared to recommendations)
+ [ ] Run a sensitivity analysis to choice of reference gene. Plan here is to
choose the gene with the least variability across samples. If the assumption
of no change in RNA composition is true, this represents a strong candidate
for a true reference gene. If the assumption is violated, this represents a
strong candidate for a "proportional gene", which is a gene whose proportion
remains the same and whose absolute amount is closely correlated to the total RNA
amount in the sample.
