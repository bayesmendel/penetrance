# Initialize Age Imputation

Initializes the age imputation process by filling missing ages with
random values between a threshold and maximum age.

## Usage

``` r
imputeAgesInit(data, threshold, max_age)
```

## Arguments

- data:

  A data frame containing family-based data

- threshold:

  Minimum age value for initialization

- max_age:

  Maximum age value for initialization

## Value

A list containing:

- data:

  The data frame with initialized ages

- na_indices:

  Indices of missing age values

## Examples

``` r
# Create sample data
data <- data.frame(
  family = c(1, 1),
  individual = c(1, 2),
  father = c(NA, 1),
  mother = c(NA, NA),
  sex = c(1, 2),
  aff = c(1, 0),
  age = c(NA, NA),
  geno = c("1/2", NA),
  isProband = c(1, 0)
)

# Initialize ages with random values between 20 and 94
result <- imputeAgesInit(data, threshold = 20, max_age = 94)

# Access the results
imputed_data <- result$data
missing_indices <- result$na_indices
```
