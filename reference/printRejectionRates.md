# Print MCMC Rejection Rates

Print MCMC Rejection Rates

## Usage

``` r
printRejectionRates(results, verbose = TRUE)
```

## Arguments

- results:

  A list of MCMC chain results.

- verbose:

  Logical, whether to print rates to console. Default is TRUE.

## Value

A named numeric vector containing the rejection rate (between 0 and 1)
for each MCMC chain. Names are of the form "Chain X" where X is the
chain number.

## Details

Extracts and prints the rejection rates from MCMC chain results.

## Examples

``` r
# Create example results list with two chains
results <- list(
  list(rejection_rate = 0.3),
  list(rejection_rate = 0.4)
)

# Get rejection rates without printing
rates <- printRejectionRates(results, verbose = FALSE)

# Print rejection rates
rates <- printRejectionRates(results)
#> Rejection rates:
#>   Chain 1: 0.30
#>   Chain 2: 0.40
```
