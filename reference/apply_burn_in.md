# Apply Burn-In

Apply Burn-In

## Usage

``` r
apply_burn_in(results, burn_in)
```

## Arguments

- results:

  A list of MCMC chain results.

- burn_in:

  The fraction roportion of results to discard as burn-in (0 to 1). The
  default is no burn-in, burn_in=0.

## Value

A list of results with burn-in applied.
