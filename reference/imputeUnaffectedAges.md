# Impute Ages for Unaffected Individuals

This function imputes ages for unaffected individuals in a dataset based
on their sex and whether they were tested, using empirical age
distributions.

## Usage

``` r
imputeUnaffectedAges(data, na_indices, empirical_density, max_age)
```

## Arguments

- data:

  A data frame containing the individual data, including columns for
  age, sex, and geno.

- na_indices:

  A vector of indices indicating the rows in the data where ages need to
  be imputed.

- empirical_density:

  A list of density objects containing the empirical density of ages for
  different groups.

- max_age:

  Integer, the maximum age considered in the analysis.

## Value

The data frame with imputed ages for unaffected individuals.
