# Make Priors

This function generates prior distributions based on user input or
default parameters. It is designed to aid in the statistical analysis of
risk proportions in populations, particularly in the context of cancer
research. The distributions are calculated for various statistical
metrics such as asymptote, threshold, median, and first quartile.

## Usage

``` r
makePriors(
  data,
  sample_size,
  ratio,
  prior_params,
  risk_proportion,
  baseline_data
)
```

## Arguments

- data:

  A data frame containing age and risk data. If NULL or contains NA
  values, default parameters are used.

- sample_size:

  Numeric, the total sample size used for risk proportion calculations.

- ratio:

  Numeric, the odds ratio (OR) or relative risk (RR) used in asymptote
  parameter calculations.

- prior_params:

  List, containing prior parameters for the beta distributions. If NULL,
  default parameters are used.

- risk_proportion:

  Data frame, with default proportions of people at risk.

- baseline_data:

  Data frame with the baseline risk data.

## Value

A list of functions representing the prior distributions for asymptote,
threshold, median, and first quartile.

## Details

The function includes internal helper functions for normalizing median
and first quartile values, and for computing beta distribution
parameters. The function handles various settings: using default
parameters, applying user inputs, and calculating parameters based on
sample size and risk proportions.

If the OR/RR ratio is provided, the asymptote parameters are computed
based on this ratio, overriding other inputs for the asymptote.

The function returns a list of distribution functions for the asymptote,
threshold, median, and first quartile, which can be used for further
statistical analysis.

## See also

[`qbeta`](https://rdrr.io/r/stats/Beta.html),
[`runif`](https://rdrr.io/r/stats/Uniform.html)
