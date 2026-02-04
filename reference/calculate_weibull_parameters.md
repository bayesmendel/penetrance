# Calculate Weibull Parameters

This function calculates the shape (`alpha`) and scale (`beta`)
parameters of a Weibull distribution given the median, first quartile,
and delta values.

## Usage

``` r
calculate_weibull_parameters(given_median, given_first_quartile, delta)
```

## Arguments

- given_median:

  The median of the data.

- given_first_quartile:

  The first quartile of the data.

- delta:

  A constant offset value.

## Value

A list containing the calculated Weibull parameters:

- alpha:

  The shape parameter of the Weibull distribution

- beta:

  The scale parameter of the Weibull distribution

## Examples

``` r
# Calculate Weibull parameters
params <- calculate_weibull_parameters(
  given_median = 50,
  given_first_quartile = 30,
  delta = 15
)
print(params)
#> $alpha
#> [1] 1.037872
#> 
#> $beta
#> [1] 49.82351
#> 
```
