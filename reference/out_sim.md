# Simulated Output Data

This dataset contains the simulated output data for the penetrance
package.

## Usage

``` r
data(out_sim)
```

## Format

A list with the following components:

- summary_stats:

  A data frame with 18000 observations of 8 variables:

  Median_Male

  :   numeric, Median value for males

  Median_Female

  :   numeric, Median value for females

  Threshold_Male

  :   numeric, Threshold value for males

  Threshold_Female

  :   numeric, Threshold value for females

  First_Quartile_Male

  :   numeric, First quartile value for males

  First_Quartile_Female

  :   numeric, First quartile value for females

  Asymptote_Male

  :   numeric, Asymptote value for males

  Asymptote_Female

  :   numeric, Asymptote value for females

- density_plots:

  A list of 1 element, mfrow: integer vector of length 2

- trace_plots:

  A list of 1 element, mfrow: integer vector of length 2

- penetrance_plot:

  A list of 2 elements: rect and text

- pdf_plots:

  A list of 2 elements: rect and text

- combined_chains:

  A list of 19 numeric vectors with 18000 elements each

- results:

  A list of 1 element which is a list of 24 elements, each with 18000
  elements

- data:

  A data frame with 4727 observations of 9 variables:

  individual

  :   integer, Individual ID

  isProband

  :   numeric, Indicator if the individual is a proband

  family

  :   integer, Family ID

  mother

  :   numeric, Mother's ID

  father

  :   numeric, Father's ID

  aff

  :   numeric, Affected status

  sex

  :   numeric, Sex of the individual

  age

  :   numeric, Age of the individual

  geno

  :   character, Genotype

## Examples

``` r
data(out_sim)
head(out_sim$summary_stats)
#>   Median_Male Median_Female Threshold_Male Threshold_Female First_Quartile_Male
#> 1    65.29099      52.24483        31.1586         14.67767            53.81554
#> 2    65.29099      52.24483        31.1586         14.67767            53.81554
#> 3    65.29099      52.24483        31.1586         14.67767            53.81554
#> 4    65.29099      52.24483        31.1586         14.67767            53.81554
#> 5    65.29099      52.24483        31.1586         14.67767            53.81554
#> 6    65.29099      52.24483        31.1586         14.67767            53.81554
#>   First_Quartile_Female Asymptote_Male Asymptote_Female
#> 1              40.01161      0.5917241        0.6169134
#> 2              40.01161      0.5917241        0.6169134
#> 3              40.01161      0.5917241        0.6169134
#> 4              40.01161      0.5917241        0.6169134
#> 5              40.01161      0.5917241        0.6169134
#> 6              40.01161      0.5917241        0.6169134
```
