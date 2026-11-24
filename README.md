
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Penetrance R package

An R package for the estimation of age-specific penetrance for complex
family-based studies in a format compatible with the Fam3PRO R package.

## Motivation

Accurate estimation of age-specific penetrance is essential for
assessing disease risk in individuals with pathogenic genetic variants
(PGVs). Penetrance refers to the probability that an individual carrying
a specific genetic variant will develop the associated trait, such as
cancer. Estimating this probability is a crucial step in clinical
decision-making and personalized risk assessment for hereditary cancer
syndromes.

The package leverages Mendelian inheritance models, which are widely
used in family-based genetic studies to assess how genetic variants are
passed down through generations. These models typically involve a
proband—an individual for whom family history and genetic data are
collected. The proband serves as the starting point for mapping out the
family’s genetic structure, including relationships and phenotypic
traits, such as cancer diagnoses. Family data, including cancer
occurrence, ages of diagnosis, and genetic test results, are collected
for the proband and their relatives. Using these data, Mendelian models
compute the likelihood of genetic configurations and disease outcomes
based on inheritance patterns.

The core methodology in the package relies on a four-parameter Weibull
distribution to model age-specific penetrance. Estimation is performed
using a Bayesian framework with Markov chain Monte Carlo (MCMC) methods,
allowing the package to provide robust and flexible penetrance
estimates. Through this approach, the package models the likelihood of
cancer occurrence across family members, even when some genotypic
information is missing or incomplete, as is common in real-world
studies.

The package also incorporates prior knowledge into the estimation
process, enabling users to specify default, custom, or study-based prior
distributions. By employing the Elston–Stewart peeling algorithm, the
package efficiently calculates likelihoods across family pedigrees,
ensuring scalability and accuracy, even in large datasets.

By providing functions for data input, prior specification, and
estimation, the package equips researchers and clinicians with tools for
estimating cancer risk in complex family-based studies. This supports
informed decision-making and preventive strategies in hereditary cancer
syndromes, where understanding the genetic basis of risk is critical for
patient care.

## Citation

If you use the `penetrance` package, please cite:

> Nicolas Kubista, Danielle Braun, Giovanni Parmigiani (2025). The
> *penetrance* R package for estimation of age specific risk in
> family-based studies. *Bioinformatics Advances*, Volume 5, Issue 1,
> vbaf154. <https://doi.org/10.1093/bioadv/vbaf154>

## Installation

To install, use

    git clone git@github.com:bayesmendel/penetrance.git

Open the source directory as a new R project and install the package
with

    devtools::install()

or install it directly from GitHub in RStudio:

    devtools::install_github("bayesmendel/penetrance")

## Quick-start guide

The following is a quick-start guide for basic usage of the package. For
greater detail on available options, please refer to the other articles.

The primary function in the package is `penetrance()`. The package
workflow includes three main parts: supplying family data and specifying
the estimation settings, estimating the posterior distribution using
MCMC, and returning samples from the approximated posterior distribution
representing the estimated penetrance function.

``` r
library(penetrance)
#> Warning: package 'penetrance' was built under R version 4.5.2
```

### Pedigree

The user must specify the `pedigree` argument as a list of data frames,
where each data frame contains one family’s data. Each data frame must
have the following columns:

- `PedigreeID`: A numeric or character identifier for the family. It
  must be consistent for all members of the family.

- `ID`: A numeric or character identifier for each individual. It must
  be unique within a pedigree.

- `Sex`: An integer where `0` indicates female and `1` indicates male.
  Unknown sex can be coded as `NA`.

- `MotherID`: The `ID` of the individual’s mother, or `NA` if the mother
  is not included in the pedigree.

- `FatherID`: The `ID` of the individual’s father, or `NA` if the father
  is not included in the pedigree.

- `isProband`: An integer where `1` indicates that the individual is a
  proband and `0` otherwise.

- `CurAge`: An integer indicating the age of censoring: current age if
  the individual is alive or age at death if deceased. Values must be
  between `1` and `max_age`; unknown ages can be coded as `NA`.

- `isAff`: An integer indicating affection status for the cancer of
  interest, with `1` for diagnosed individuals and `0` for unaffected
  individuals. Unknown status can be coded as `NA`.

- `Age`: An integer indicating age at cancer diagnosis. It should be
  `NA` if `isAff` is `0` or `NA`; otherwise it must be between `1` and
  `max_age` and should not exceed `CurAge`.

- `Geno`: An integer indicating germline genetic test results, with `1`
  for carriers and `0` for non-carriers. Unknown or untested individuals
  can be coded as `NA`.

### Model specification

Available options include:

- `pedigree`: A list of data frames containing pedigrees in the format
  described above.
- `twins`: A list specifying identical twins or triplets. For example,
  `list(c("ora024", "ora027"), c("aey063", "aey064"))` specifies two
  identical-twin pairs.
- `n_chains`: Integer, the number of chains for parallel computation.
  Default is 1.
- `n_iter_per_chain`: Integer, the number of iterations for each chain.
  Default is 10000.
- `ncores`: Integer, the number of cores for parallel computation.
  Default is 6.
- `baseline_data`: Absolute age-specific baseline risk data. With
  `sex_specific = TRUE`, this is a data frame with `Male` and `Female`
  columns; with `sex_specific = FALSE`, it is a numeric vector or
  single-column data frame.
- `max_age`: Integer, the maximum age considered for analysis. Default
  is 94.
- `remove_proband`: Logical, indicating whether to remove probands from
  the analysis. Default is `FALSE`.
- `age_imputation`: Logical, indicating whether to perform age
  imputation. Default is `FALSE`.
- `median_max`: Logical, indicating whether to use the baseline median
  age or `max_age` as an upper bound for the median proposal. Default is
  `TRUE`.
- `BaselineNC`: Logical, indicating that non-carrier penetrance is
  assumed to equal baseline penetrance. Only `TRUE` is currently
  supported.
- `var`: Numeric vector containing proposal variances for the
  Metropolis–Hastings algorithm.
- `burn_in`: Numeric, the fraction of results to discard as burn-in.
  Default is 0.
- `thinning_factor`: Integer, the factor by which to thin the results.
  Default is 1.
- `imp_interval`: Integer, the interval at which age imputation is
  performed when `age_imputation = TRUE`. Default is 100.
- `distribution_data`: Data used to generate prior distributions.
- `allele_freq`: Numeric, the population allele frequency of the risk
  variant. Default is 0.0001.
- `sample_size`: Optional numeric sample size for distribution
  generation.
- `ratio`: Optional numeric ratio for distribution generation.
- `prior_params`: A list containing prior-distribution parameters.
- `risk_proportion`: Numeric, the risk proportion used for distribution
  generation.
- `summary_stats`: Logical, indicating whether to include summary
  statistics in the output. Default is `TRUE`.
- `rejection_rates`: Logical, indicating whether to include rejection
  rates in the output. Default is `TRUE`.
- `density_plots`: Logical, indicating whether to include density plots
  in the output. Default is `TRUE`.
- `plot_trace`: Logical, indicating whether to include trace plots in
  the output. Default is `TRUE`.
- `penetrance_plot`: Logical, indicating whether to include cumulative
  penetrance plots in the output. Default is `TRUE`.
- `penetrance_plot_pdf`: Logical, indicating whether to include
  penetrance density plots in the output. Default is `TRUE`.
- `plot_loglikelihood`: Logical, indicating whether to include
  log-likelihood plots in the output. Default is `TRUE`.
- `plot_acf`: Logical, indicating whether to include autocorrelation
  plots for posterior samples. Default is `TRUE`.
- `probCI`: Numeric, the probability level for credible intervals in
  penetrance plots. Default is 0.95.
- `sex_specific`: Logical, indicating whether to use sex-specific
  parameters in the analysis. Default is `TRUE`.

### Additional user inputs

- The `penetrance()` function takes absolute age-specific probabilities
  of developing cancer as the `baseline_data` input. With the default
  `BaselineNC = TRUE`, this baseline is assumed to represent non-carrier
  penetrance. For rare variants, this is considered a reasonable
  assumption.

- The allele frequency (`allele_freq`) defaults to 0.0001. The function
  converts it to carrier prevalence using the approximation
  `carrier_prevalence` ≈ 2 × `allele_freq` for rare autosomal dominant
  conditions.

- The `penetrance()` function includes an option for automatic age
  imputation through `age_imputation`. Ages are imputed during the MCMC
  routine based on affection status, sex, and degree of relationship to
  the proband who carries the pathogenic variant.

- Monozygotic twins or triplets can be specified using the `twins`
  argument. For example:

``` r
twins <- list(c("ora024", "ora027"), c("aey063", "aey064"))
```
