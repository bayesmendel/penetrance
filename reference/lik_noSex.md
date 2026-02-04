# Likelihood Calculation without Sex Differentiation

This function calculates the likelihood for an individual based on
Weibull distribution parameters without considering sex differentiation.

## Usage

``` r
lik_noSex(
  i,
  data,
  alpha,
  beta,
  delta,
  gamma,
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
  Must include columns for 'age', 'aff' (affection status), and 'geno'
  (genotype).

- alpha:

  Numeric, Weibull distribution shape parameter.

- beta:

  Numeric, Weibull distribution scale parameter.

- delta:

  Numeric, shift parameter for the Weibull function.

- gamma:

  Numeric, asymptote parameter (only scales the entire distribution).

- max_age:

  Integer, maximum age considered in the analysis.

- baselineRisk:

  Numeric vector, baseline risk for each age.

- BaselineNC:

  Logical, indicates if non-carrier penetrance should be based on SEER
  data or the calculated non-carrier penetrance.

- prev:

  Numeric, prevalence of the risk allele in the population.

## Value

Numeric vector, containing likelihood values for unaffected and affected
individuals.
