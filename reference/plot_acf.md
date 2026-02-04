# Plot Autocorrelation for Multiple MCMC Chains (Posterior Samples)

This function plots the autocorrelation for sex-specific or
non-sex-specific posterior samples across multiple MCMC chains. It
defaults to key parameters like `asymptote_male_samples`,
`asymptote_female_samples`, etc.

## Usage

``` r
plot_acf(results, n_chains, max_lag = 50)
```

## Arguments

- results:

  A list of MCMC chain results.

- n_chains:

  The number of chains.

- max_lag:

  Integer, the maximum lag to be considered for the autocorrelation
  plot. Default is 50.

## Value

A series of autocorrelation plots for each chain.
