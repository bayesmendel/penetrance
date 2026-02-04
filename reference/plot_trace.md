# Plot MCMC Trace Plots

Plot MCMC Trace Plots

## Usage

``` r
plot_trace(results, n_chains, verbose = FALSE)
```

## Arguments

- results:

  A list of MCMC chain results.

- n_chains:

  Integer, the number of chains.

- verbose:

  Logical, whether to print progress messages. Default is FALSE.

## Value

No return value, called for side effects. Creates trace plots for each
parameter.

## Examples

``` r
# Create example results list
results <- list(
  list(
    median_samples = rnorm(100),
    threshold_samples = runif(100),
    first_quartile_samples = rgamma(100, 2, 2),
    asymptote_samples = rbeta(100, 2, 2)
  )
)

# Generate trace plots
plot_trace(results, n_chains = 1)
```
