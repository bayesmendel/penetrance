# Calculate Empirical Age Density

Calculates empirical age density distributions for different subgroups
in the data, separated by sex and genetic testing status.

## Usage

``` r
calculateEmpiricalDensity(
  data,
  aff_column = "aff",
  age_column = "age",
  sex_column = "sex",
  geno_column = "geno",
  n_points = 10000,
  sex_specific = TRUE
)
```

## Arguments

- data:

  A data frame containing the family data

- aff_column:

  Name of the affection status column

- age_column:

  Name of the age column

- sex_column:

  Name of the sex column

- geno_column:

  Name of the genotype column

- n_points:

  Number of points to use in density estimation

- sex_specific:

  Logical; whether to calculate sex-specific densities

## Value

A list of density objects for different subgroups (tested/untested,
male/female)
