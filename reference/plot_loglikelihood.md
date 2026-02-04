# Plot Log-Likelihood for Multiple MCMC Chains

This function plots the log-likelihood values across iterations for
multiple MCMC chains. It helps visualize the convergence of the chains
based on the log-likelihood values.

## Usage

``` r
plot_loglikelihood(results, n_chains)
```

## Arguments

- results:

  A list of MCMC chain results, each containing the
  `loglikelihood_current` values.

- n_chains:

  The number of chains.

## Value

A series of log-likelihood plots for each chain.
