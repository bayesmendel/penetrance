# Calculate Log Likelihood without Sex Differentiation

This function calculates the log likelihood for a set of parameters and
data without considering sex differentiation using the clipp package.

## Usage

``` r
mhLogLikelihood_clipp_noSex(
  paras,
  families,
  twins,
  max_age,
  baseline_data,
  prev,
  geno_freq,
  trans,
  BaselineNC,
  ncores
)
```

## Arguments

- paras:

  Numeric vector, the parameters for the Weibull distribution and
  scaling factors. Should contain in order: gamma, delta, given_median,
  given_first_quartile.

- families:

  Data frame, containing pedigree information with columns for 'age',
  'aff' (affection status), and 'geno' (genotype).

- twins:

  Information on monozygous twins or triplets in the pedigrees.

- max_age:

  Integer, maximum age considered in the analysis.

- baseline_data:

  Numeric vector, baseline risk data for each age.

- prev:

  Numeric, prevalence of the risk allele in the population.

- geno_freq:

  Numeric vector, represents the frequency of the risk type and its
  complement in the population.

- trans:

  Numeric matrix, transition matrix that defines the probabilities of
  allele transmission from parents to offspring.

- BaselineNC:

  Logical, indicates if non-carrier penetrance should be based on the
  baseline data or the calculated non-carrier penetrance.

- ncores:

  Integer, number of cores to use for parallel computation.

## Value

Numeric, the calculated log likelihood.

## References

Details about the clipp package and methods can be found in the package
documentation.
