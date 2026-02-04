# Combine Chains for Non-Sex-Specific Estimation

Combines the posterior samples from multiple MCMC chains for
non-sex-specific estimations.

## Usage

``` r
combine_chains_noSex(results)
```

## Arguments

- results:

  A list of MCMC chain results, where each element contains posterior
  samples of parameters.

## Value

A list with combined results, including samples for median, threshold,
first quartile, asymptote values, log-likelihoods, and log-acceptance
ratios.
