# penetrance: A Package for Penetrance Estimation

A comprehensive package for penetrance estimation in family-based
studies. This package implements Bayesian methods using
Metropolis-Hastings algorithm for estimating age-specific penetrance of
genetic variants. It supports both sex-specific and non-sex-specific
analyses, and provides various visualization tools for examining MCMC
results.

This function implements the Independent Metropolis-Hastings algorithm
for Bayesian penetrance estimation of cancer risk. It utilizes parallel
computing to run multiple chains and provides various options for
analyzing and visualizing the results.

## Usage

``` r
penetrance(
  pedigree,
  twins = NULL,
  n_chains = 1,
  n_iter_per_chain = 10000,
  ncores = 6,
  max_age = 94,
  baseline_data = baseline_data_default,
  remove_proband = FALSE,
  age_imputation = FALSE,
  median_max = TRUE,
  BaselineNC = TRUE,
  var = c(0.1, 0.1, 2, 2, 5, 5, 5, 5),
  burn_in = 0,
  thinning_factor = 1,
  imp_interval = 100,
  distribution_data = distribution_data_default,
  prev = 1e-04,
  sample_size = NULL,
  ratio = NULL,
  prior_params = prior_params_default,
  risk_proportion = risk_proportion_default,
  summary_stats = TRUE,
  rejection_rates = TRUE,
  density_plots = TRUE,
  plot_trace = TRUE,
  penetrance_plot = TRUE,
  penetrance_plot_pdf = TRUE,
  plot_loglikelihood = TRUE,
  plot_acf = TRUE,
  probCI = 0.95,
  sex_specific = TRUE
)
```

## Arguments

- pedigree:

  A list of data frames, where each data frame represents a single
  pedigree and contains the following columns:

  - `PedigreeID`: A numeric or character identifier for the
    family/pedigree. Must be consistent for all members of the same
    family within a data frame.

  - `ID`: A unique numeric or character identifier for each individual
    within their respective pedigree data frame.

  - `Sex`: An integer representing biological sex: `0` for female, `1`
    for male. Use `NA` for unknown sex.

  - `MotherID`: The `ID` of the individual's mother. Should correspond
    to an `ID` within the same pedigree data frame or be `NA` if the
    mother is not in the pedigree (founder).

  - `FatherID`: The `ID` of the individual's father. Should correspond
    to an `ID` within the same pedigree data frame or be `NA` if the
    father is not in the pedigree (founder).

  - `isProband`: An integer indicating if the individual is a proband:
    `1` for proband, `0` otherwise.

  - `CurAge`: An integer representing the age of censoring. This is the
    current age if the individual is alive, or the age at death if
    deceased. Must be between `1` and `max_age`. Use `NA` for unknown
    ages (but note this may affect analysis or require imputation).

  - `isAff`: An integer indicating the affection status for the cancer
    of interest: `1` if diagnosed, `0` if unaffected. Use `NA` for
    unknown status.

  - `Age`: An integer representing the age at cancer diagnosis. Should
    be `NA` if `isAff` is `0` or `NA`. Must be between `1` and
    `max_age`, and less than or equal to `CurAge`. Use `NA` for unknown
    diagnosis age (but note this may affect analysis or require
    imputation).

  - `Geno`: An integer representing the germline genetic test result:
    `1` for carrier (positive), `0` for non-carrier (negative). Use `NA`
    for unknown or untested individuals.

- twins:

  A list specifying identical twins or triplets in the family. Each
  element of the list should be a vector containing the `ID`s of the
  identical siblings within a pedigree. For example:
  `list(c("ID1", "ID2"), c("ID3", "ID4", "ID5"))`. Default is `NULL`.

- n_chains:

  Integer, the number of chains for parallel computation. Default is 1.

- n_iter_per_chain:

  Integer, the number of iterations for each chain. Default is 10000.

- ncores:

  Integer, the number of cores for parallel computation. Default is 6.

- max_age:

  Integer, the maximum age considered for analysis. Default is 94.

- baseline_data:

  Data providing the absolute age-specific baseline risk (probability)
  of developing the cancer in the general population (e.g., from SEER
  database). All probability values must be between 0 and 1. - If
  `sex_specific = TRUE` (default): A data frame with columns 'Male' and
  'Female', where each column contains the age-specific probabilities
  for that sex. The number of rows should ideally correspond to
  `max_age`. - If `sex_specific = FALSE`: A numeric vector or a
  single-column data frame containing the age-specific probabilities for
  the combined population. The length (or number of rows) should ideally
  correspond to `max_age`. Default data is provided for Colorectal
  cancer from SEER (up to age 94). If the number of rows/length does not
  match `max_age`, the data will be truncated or extended with the last
  value.

- remove_proband:

  Logical, indicating whether to remove probands from the analysis.
  Default is FALSE.

- age_imputation:

  Logical, indicating whether to perform age imputation. Default is
  FALSE.

- median_max:

  Logical, indicating whether to use the baseline median age or
  `max_age` as an upper bound for the median proposal. Default is TRUE.

- BaselineNC:

  Logical, indicating that the non-carrier penetrance is assumed to be
  the baseline penetrance. Default is TRUE.

- var:

  Numeric vector, variances for the proposal distribution in the
  Metropolis-Hastings algorithm. Default is
  `c(0.1, 0.1, 2, 2, 5, 5, 5, 5)`.

- burn_in:

  Numeric, the fraction of results to discard as burn-in (0 to 1).
  Default is 0 (no burn-in).

- thinning_factor:

  Integer, the factor by which to thin the results. Default is 1 (no
  thinning).

- imp_interval:

  Integer, the interval at which age imputation should be performed when
  age_imputation = TRUE.

- distribution_data:

  Data for generating prior distributions.

- prev:

  Numeric, prevalence of the carrier status. Default is 0.0001.

- sample_size:

  Optional numeric, sample size for distribution generation.

- ratio:

  Optional numeric, ratio parameter for distribution generation.

- prior_params:

  List, parameters for prior distributions.

- risk_proportion:

  Numeric, proportion of risk for distribution generation.

- summary_stats:

  Logical, indicating whether to include summary statistics in the
  output. Default is TRUE.

- rejection_rates:

  Logical, indicating whether to include rejection rates in the output.
  Default is TRUE.

- density_plots:

  Logical, indicating whether to include density plots in the output.
  Default is TRUE.

- plot_trace:

  Logical, indicating whether to include trace plots in the output.
  Default is TRUE.

- penetrance_plot:

  Logical, indicating whether to include penetrance plots in the output.
  Default is TRUE.

- penetrance_plot_pdf:

  Logical, indicating whether to include PDF plots in the output.
  Default is TRUE.

- plot_loglikelihood:

  Logical, indicating whether to include log-likelihood plots in the
  output. Default is TRUE.

- plot_acf:

  Logical, indicating whether to include autocorrelation function (ACF)
  plots for posterior samples. Default is TRUE.

- probCI:

  Numeric, probability level for credible intervals in penetrance plots.
  Must be between 0 and 1. Default is 0.95.

- sex_specific:

  Logical, indicating whether to use sex-specific parameters in the
  analysis. Default is TRUE.

## Value

A list containing combined results from all chains, including optional
statistics and plots.

## Details

Key features:

- Bayesian estimation of penetrance using family-based data

- Support for sex-specific and non-sex-specific analyses

- Age imputation for missing data

- Visualization tools for MCMC diagnostics

- Integration with the clipp package for likelihood calculations

## See also

Useful links:

- <https://github.com/nicokubi/penetrance>

- <https://nicokubi.github.io/penetrance/>

## Author

**Maintainer**: Sol Rosito <bmendel@jimmy.harvard.edu>

Authors:

- BayesMendel Lab

## Examples

``` r
# Create example baseline data (simplified for demonstration)
baseline_data_default <- data.frame(
  Age = 1:94,
  Female = rep(0.01, 94),
  Male = rep(0.01, 94)
)

# Create example distribution data
distribution_data_default <- data.frame(
  Age = 1:94,
  Risk = rep(0.01, 94)
)

# Create example prior parameters
prior_params_default <- list(
  shape = 2,
  scale = 50
)

# Create example risk proportion
risk_proportion_default <- 0.5

# Create a simple example pedigree
example_pedigree <- data.frame(
  PedigreeID = rep(1, 4),
  ID = 1:4,
  Sex = c(1, 0, 1, 0),  # 1 for male, 0 for female
  MotherID = c(NA, NA, 2, 2),
  FatherID = c(NA, NA, 1, 1),
  isProband = c(0, 0, 1, 0),
  CurAge = c(70, 68, 45, 42),
  isAff = c(0, 0, 1, 0),
  Age = c(NA, NA, 40, NA),
  Geno = c(NA, NA, 1, NA)
)

# Basic usage with minimal iterations
result <- penetrance(
  pedigree = list(example_pedigree),
  n_chains = 1,
  n_iter_per_chain = 10,  # Very small number for example
  ncores = 1,             # Single core for example
  summary_stats = TRUE,
  plot_trace = FALSE,     # Disable plots for quick example
  density_plots = FALSE,
  penetrance_plot = FALSE,
  penetrance_plot_pdf = FALSE,
  plot_loglikelihood = FALSE,
  plot_acf = FALSE
)
#> Rejection rates:
#>   Chain 1: 0.20

# View basic results
head(result$summary_stats)
#>   Median_Male    Median_Female   Threshold_Male  Threshold_Female
#>  Min.   :40.00   Min.   :50.00   Min.   :32.99   Min.   :25.00   
#>  1st Qu.:40.66   1st Qu.:50.61   1st Qu.:33.44   1st Qu.:25.95   
#>  Median :42.86   Median :52.98   Median :33.63   Median :26.32   
#>  Mean   :42.33   Mean   :52.38   Mean   :33.69   Mean   :26.27   
#>  3rd Qu.:43.34   3rd Qu.:53.20   3rd Qu.:33.81   3rd Qu.:26.60   
#>  Max.   :44.74   Max.   :54.92   Max.   :35.00   Max.   :27.21   
#>  First_Quartile_Male First_Quartile_Female Asymptote_Male    Asymptote_Female
#>  Min.   :38.01       Min.   :31.07         Min.   :0.05812   Min.   :0.5086  
#>  1st Qu.:39.12       1st Qu.:32.52         1st Qu.:0.11858   1st Qu.:0.5667  
#>  Median :40.03       Median :32.94         Median :0.13965   Median :0.5962  
#>  Mean   :39.61       Mean   :32.60         Mean   :0.14527   Mean   :0.6219  
#>  3rd Qu.:40.14       3rd Qu.:33.16         3rd Qu.:0.16854   3rd Qu.:0.7064  
#>  Max.   :40.66       Max.   :33.48         Max.   :0.23874   Max.   :0.7546  
```
