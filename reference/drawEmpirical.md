# Draw Ages Using the Inverse CDF Method from Empirical Density

This function draws ages using the inverse CDF method from empirical
density data, based on sex and whether the individual was tested.

## Usage

``` r
drawEmpirical(empirical_density, sex, tested, sex_specific = TRUE)
```

## Arguments

- empirical_density:

  A list of density objects containing the empirical density of ages for
  different groups.

- sex:

  Numeric, the sex of the individual (1 for male, 2 for female).

- tested:

  Logical, indicating whether the individual was tested (has a non-NA
  'geno' value).

- sex_specific:

  Logical, indicating whether the imputation should be sex-specific.
  Default is TRUE.

## Value

A single age value drawn from the appropriate empirical density data.
