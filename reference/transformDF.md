# Transform Data Frame

This function transforms a data frame from the standard format used in
PanelPRO into the required format which conforms to the requirements of
penetrance (and clipp).

## Usage

``` r
transformDF(df)
```

## Arguments

- df:

  The input data frame in the usual PanelPRO format.

## Value

A data frame in the format required for clipp with the following
columns:

- individual:

  ID of the individual

- isProband:

  Indicator if the individual is a proband

- family:

  Family ID

- mother:

  Mother's ID

- father:

  Father's ID

- aff:

  Affection status

- sex:

  Sex (2 for female, 1 for male)

- age:

  Age at diagnosis or current age

- geno:

  Genotype information

## Examples

``` r
# Create example data frame
df <- data.frame(
  ID = 1:2,
  PedigreeID = c(1,1),
  Sex = c(0,1),
  MotherID = c(NA,1),
  FatherID = c(NA,NA),
  isProband = c(1,0),
  CurAge = c(45,20),
  isAff = c(1,0),
  Age = c(40,NA),
  Geno = c(1,0)
)

# Transform the data frame
transformed_df <- transformDF(df)
```
