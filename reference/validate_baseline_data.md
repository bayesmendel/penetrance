# Validate Baseline Penetrance Data

This function validates baseline penetrance data to ensure it represents
age-specific probabilities rather than cumulative risk. It checks for
monotonicity and whether the sum exceeds 1, both of which suggest the
data may be cumulative rather than age-specific.

## Usage

``` r
validate_baseline_data(
  baseline_data,
  sex_specific = TRUE,
  param_name = "baseline_data",
  tolerance = 1e-10
)
```

## Arguments

- baseline_data:

  The baseline data to validate. Can be:

  - A data frame with 'Male' and 'Female' columns (when sex_specific =
    TRUE)

  - A numeric vector (when sex_specific = FALSE)

  - A single-column data frame (when sex_specific = FALSE)

- sex_specific:

  Logical, indicating whether the data is sex-specific. Default is TRUE.

- param_name:

  Character string specifying the parameter name (for messages). Default
  is "baseline_data".

- tolerance:

  Numeric value for checking strict monotonicity (to account for
  floating point precision). Default is 1e-10.

## Value

Logical value TRUE if validation passes (with possible warnings),
otherwise stops with an error.

## Details

The function performs the following checks:

- For monotonicity: If values are strictly non-decreasing (monotonically
  increasing), this suggests cumulative risk rather than age-specific
  probabilities. A warning is issued.

- For sum \> 1: If the sum of all probabilities exceeds 1, this is
  problematic because these should be age-specific probabilities. A
  warning is issued.

- Individual values must be between 0 and 1 (probabilities)

- No NA or infinite values are allowed

Age-specific baseline risk represents the probability of developing
disease at each specific age, while cumulative risk represents the total
probability up to that age. For proper penetrance estimation,
age-specific (not cumulative) risk should be used.

## Examples

``` r
# Valid age-specific data (varies, not monotone)
age_specific <- c(0.001, 0.002, 0.003, 0.002, 0.004, 0.003, 0.005)
validate_baseline_data(age_specific, sex_specific = FALSE)
#> [1] TRUE

# Valid sex-specific data
baseline_df <- data.frame(
  Male = c(0.001, 0.002, 0.001, 0.003),
  Female = c(0.002, 0.003, 0.002, 0.004)
)
validate_baseline_data(baseline_df, sex_specific = TRUE)
#> [1] TRUE

if (FALSE) { # \dontrun{
# Will trigger warnings
# Monotone increasing (suggests cumulative risk)
cumulative <- c(0.001, 0.002, 0.003, 0.004, 0.005)
validate_baseline_data(cumulative, sex_specific = FALSE)

# Sum greater than 1
high_values <- rep(0.1, 15)  # sum = 1.5
validate_baseline_data(high_values, sex_specific = FALSE)

# Invalid data
invalid_data <- c(0.001, -0.002, 0.003)  # Negative value
validate_baseline_data(invalid_data, sex_specific = FALSE)
} # }
```
