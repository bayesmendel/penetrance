# Generate Summary for Non-Sex-Specific Estimation

Generates summary statistics for the combined MCMC results for
non-sex-specific estimations.

## Usage

``` r
generate_summary_noSex(data, verbose = FALSE)
```

## Arguments

- data:

  A list containing combined results of MCMC chains, typically the
  output of `combine_chains_noSex`.

- verbose:

  Logical, whether to print summary to console. Default is FALSE.

## Value

A data.frame containing summary statistics (min, 1st quartile, median,
mean, 3rd quartile, max) for median, threshold, first quartile, and
asymptote values.
