# Validate Weibull Parameters

This function validates the given parameters for calculating Weibull
distribution.

## Usage

``` r
validate_weibull_parameters(
  given_first_quartile,
  given_median,
  threshold,
  asymptote
)
```

## Arguments

- given_first_quartile:

  The first quartile of the data.

- given_median:

  The median of the data.

- threshold:

  A constant threshold value.

- asymptote:

  A constant asymptote value (gamma).

## Value

Logical value indicating whether the parameters are valid (TRUE) or not
(FALSE)

## Examples

``` r
# Validate parameters
is_valid <- validate_weibull_parameters(
  given_first_quartile = 30,
  given_median = 50,
  threshold = 15,
  asymptote = 0.8
)
print(is_valid)
#> [1] TRUE
```
