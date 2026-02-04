# Plot Weibull Probability Density Function with Credible Intervals

This function plots the Weibull PDF with credible intervals for the
given MCMC results. It allows for visualization of density curves based
on the posterior samples.

## Usage

``` r
plot_pdf(combined_chains, prob, max_age, sex = "NA")
```

## Arguments

- combined_chains:

  List of combined MCMC chain results containing posterior samples for
  penetrance parameters.

- prob:

  Numeric, probability level for confidence intervals (between 0 and 1).

- max_age:

  Integer, maximum age to plot.

- sex:

  Character, one of "Male", "Female", or "NA" for non-sex-specific.
  Default is "NA".

## Value

A plot showing the Weibull PDF with credible intervals.
