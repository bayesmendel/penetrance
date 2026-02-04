# Penetrance Function

Calculates the penetrance for an individual based on Weibull
distribution parameters. This function estimates the probability of
developing cancer given the individual's genetic and demographic
information.

## Usage

``` r
lik.fn(
  i,
  data,
  alpha_male,
  alpha_female,
  beta_male,
  beta_female,
  delta_male,
  delta_female,
  gamma_male,
  gamma_female,
  max_age,
  baselineRisk,
  BaselineNC,
  prev
)
```

## Arguments

- i:

  Integer, index of the individual in the data set.

- data:

  Data frame, containing individual demographic and genetic information.
  Must include columns for 'sex', 'age', 'aff' (affection status), and
  'geno' (genotype).

- alpha_male:

  Numeric, Weibull distribution shape parameter for males.

- alpha_female:

  Numeric, Weibull distribution shape parameter for females.

- beta_male:

  Numeric, Weibull distribution scale parameter for males.

- beta_female:

  Numeric, Weibull distribution scale parameter for females.

- delta_male:

  Numeric, shift parameter for the Weibull function for males.

- delta_female:

  Numeric, shift parameter for the Weibull function for females.

- gamma_male:

  Numeric, asymptote parameter for males (only scales the entire
  distribution).

- gamma_female:

  Numeric, asymptote parameter for females (only scales the entire
  distribution).

- max_age:

  Integer, maximum age considered in the analysis.

- baselineRisk:

  Numeric matrix, baseline risk for each age by sex. Rows correspond to
  sex (1 for male, 2 for female) and columns to age.

- BaselineNC:

  Logical, indicates if non-carrier penetrance should be based on SEER
  data.

- prev:

  Numeric, prevalence of the risk allele in the population.

## Value

Numeric vector, containing penetrance values for unaffected and affected
individuals.
