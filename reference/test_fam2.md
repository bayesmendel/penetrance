# Processed Family Data

A dataset containing processed information about the first simulated 130
families. These families are referenced in the vigniette
simulation_study_real.Rmd The user must specify the `pedigree` argument
as a data frame which contains the family data (see `test_fam`). The
family data must be in the correct format with the following columns:

## Usage

``` r
test_fam2
```

## Format

A list of processed family data.

## Source

Generated for package example

## Details

- `ID`:

  A numeric value representing the unique identifier for each
  individual. There should be no duplicated entries.

- `Sex`:

  A numeric value where `0` indicates female and `1` indicates male.
  Missing entries are not currently supported.

- `MotherID`:

  A numeric value representing the unique identifier for an individual's
  mother.

- `FatherID`:

  A numeric value representing the unique identifier for an individual's
  father.

- `isProband`:

  A numeric value where `1` indicates the individual is a proband and
  `0` otherwise.

- `CurAge`:

  A numeric value indicating the age of censoring (current age if the
  person is alive or age at death if the person is deceased). Allowed
  ages range from `1` to `94`.

- `isAff`:

  A numeric value indicating the affection status of cancer, with `1`
  for diagnosed individuals and `0` otherwise. Missing entries are not
  supported.

- `Age`:

  A numeric value indicating the age of cancer diagnosis, encoded as
  `NA` if the individual was not diagnosed. Allowed ages range from `1`
  to `94`.

- `Geno`:

  A column for germline testing or tumor marker testing results.
  Positive results should be coded as `1`, negative results as `0`, and
  unknown results as `NA` or left empty.
