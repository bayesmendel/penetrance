# Generate Posterior Density Plots

Generates histograms of the posterior samples for the different
parameters

## Usage

``` r
generate_density_plots(data)
```

## Arguments

- data:

  A list with combined results.

## Value

No return value, called for side effects. Creates density plots for each
parameter.

## Examples

``` r
# Create example data
data <- list(
  median_male_results = rnorm(1000, 50, 5),
  median_female_results = rnorm(1000, 45, 5),
  threshold_male_results = runif(1000, 20, 30),
  threshold_female_results = runif(1000, 25, 35),
  asymptote_male_results = rbeta(1000, 2, 2),
  asymptote_female_results = rbeta(1000, 2, 2)
)

# Generate density plots
old_par <- par(no.readonly = TRUE)  # Save old par settings
generate_density_plots(data)

par(old_par)  # Restore old par settings
```
