# Calculate Baseline Risk

This function extracts the penetrance data for a specified cancer type,
gene, race, and penetrance type from the provided database.

## Usage

``` r
calculateBaseline(cancer_type, gene, race, type, db)
```

## Arguments

- cancer_type:

  The type of cancer for which the risk is being calculated.

- gene:

  The gene of interest for which the risk is being calculated.

- race:

  The race of the individual.

- type:

  The type of penetrance calculation.

- db:

  The dataset used for the calculation, containing penetrance data.

## Value

A matrix of penetrance data for the specified parameters.
