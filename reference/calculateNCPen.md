# Calculate Age-Specific Non-Carrier Penetrance

This function calculates the age-specific non-carrier penetrance based
on SEER baseline data, penetrances for carriers, and allele frequencies.
It adjusts penetrance estimates for genetic testing by incorporating the
genetic risk attributable to specified alleles.

## Usage

``` r
calculateNCPen(SEER_baseline, alpha, beta, delta, gamma, prev, max_age)
```

## Arguments

- SEER_baseline:

  Numeric, the baseline penetrance derived from SEER data for the
  general population without considering genetic risk factors.

- alpha:

  Numeric, shape parameter for the Weibull distribution used to model
  carrier risk.

- beta:

  Numeric, scale parameter for the Weibull distribution used to model
  carrier risk.

- delta:

  Numeric, location parameter for the Weibull distribution used to model
  carrier risk.

- gamma:

  Numeric, scaling factor applied to the Weibull distribution to adjust
  carrier risk.

- prev:

  Numeric, the prevalence of the risk allele in the population.

- max_age:

  Integer, the maximum age up to which the calculations are performed.

## Value

A list containing:

- weightedCarrierRisk:

  Numeric vector, the weighted risk for carriers at each age based on
  prevalence.

- yearlyProb:

  Numeric vector, the yearly probability of not getting the disease at
  each age.

- cumulativeProb:

  Numeric vector, the cumulative probability of not getting the disease
  up to each age.
