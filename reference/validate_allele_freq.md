# Validate Allele Frequency Input

This function validates whether the provided input is a valid allele
frequency. Allele frequencies must be numeric, scalar, and within the
range 0, 1. The function also provides warnings for edge cases and
common mistakes.

## Usage

``` r
validate_allele_freq(
  allele_freq,
  param_name = "allele_freq",
  warn_threshold = 0.01
)
```

## Arguments

- allele_freq:

  The allele frequency value to validate. Should be a numeric value
  between 0 and 1.

- param_name:

  Character string specifying the parameter name (for error messages).
  Default is "allele_freq".

- warn_threshold:

  Numeric value above which to issue a warning about unusually high
  allele frequency. Default is 0.01 (1%).

## Value

Logical value TRUE if the allele frequency is valid, otherwise stops
with an error message.

## Details

The function checks:

- Whether the input is numeric

- Whether the input is a single value (not a vector)

- Whether the value is between 0 and 1 (inclusive)

- Whether the value is unusually high (\> warn_threshold), which may
  indicate the user provided carrier prevalence instead of allele
  frequency

- Whether the value is exactly 0 or 1, which may not be biologically
  meaningful

## Examples

``` r
# Valid allele frequencies
validate_allele_freq(0.0001)  # Common for rare variants
#> [1] TRUE
validate_allele_freq(0.001)
#> [1] TRUE
validate_allele_freq(0.05)
#> Warning: Warning: 'allele_freq' is 0.05 (5%), which is relatively high for a disease-associated variant. Please verify that this is the allele frequency (p) and not the carrier prevalence (approximately 2p). For example, if the carrier prevalence is 2%, the allele frequency should be approximately 1% (0.01).
#> [1] TRUE

if (FALSE) { # \dontrun{
# Invalid inputs (will throw errors)
validate_allele_freq("0.001")       # Not numeric
validate_allele_freq(c(0.001, 0.002))  # Vector instead of scalar
validate_allele_freq(-0.001)       # Negative value
validate_allele_freq(1.5)          # Greater than 1
validate_allele_freq(NA)           # Missing value

# Valid but will trigger warnings
validate_allele_freq(0.02)         # Unusually high (>1%), warning issued
validate_allele_freq(0)            # Edge case, warning issued
validate_allele_freq(1)            # Edge case, warning issued
} # }
```
