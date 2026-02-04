# Apply Thinning

Apply Thinning

## Usage

``` r
apply_thinning(results, thinning_factor)
```

## Arguments

- results:

  A list of MCMC chain results.

- thinning_factor:

  The factor by which to thin the results (positive integer). The
  default thinning factor is 1, which implies no thinning.

## Value

A list of results with thinning applied.
