
<!-- README.md is generated from README.Rmd. Please edit that file -->

# genesengStats <a href=#><img src='inst/logo.png' align="right" height="139" /></a>

> Basic statistical methods to make preliminary analysis

`genesengStats` is a dependency of Geneseng app used to quickly describe
tabular datasets.

## Installation

The latest version can be installed from GitHub as follows:

``` r
install.packages("devtools")
devtools::install_github("geneseng/genesengStats")
```

## Example

``` r
suppressWarnings({
  genesengStats::geneseng_summary_stats(data = iris, group = "Species")[[1]]
})
#>       biomarker      group  n n distinct min median mean   sd  iqr max NA's
#> 1  Petal.Length     setosa 50          9 1.0   1.50 1.46 0.17 0.18 1.9    0
#> 2  Petal.Length versicolor 50         19 3.0   4.35 4.26 0.47 0.60 5.1    0
#> 3  Petal.Length  virginica 50         20 4.5   5.55 5.55 0.55 0.78 6.9    0
#> 4   Petal.Width     setosa 50          6 0.1   0.20 0.25 0.11 0.10 0.6    0
#> 5   Petal.Width versicolor 50          9 1.0   1.30 1.33 0.20 0.30 1.8    0
#> 6   Petal.Width  virginica 50         12 1.4   2.00 2.03 0.27 0.50 2.5    0
#> 7  Sepal.Length     setosa 50         15 4.3   5.00 5.01 0.35 0.40 5.8    0
#> 8  Sepal.Length versicolor 50         21 4.9   5.90 5.94 0.52 0.70 7.0    0
#> 9  Sepal.Length  virginica 50         21 4.9   6.50 6.59 0.64 0.67 7.9    0
#> 10  Sepal.Width     setosa 50         16 2.3   3.40 3.43 0.38 0.48 4.4    0
#> 11  Sepal.Width versicolor 50         14 2.0   2.80 2.77 0.31 0.48 3.4    0
#> 12  Sepal.Width  virginica 50         13 2.2   3.00 2.97 0.32 0.38 3.8    0
#>    Shapiro's test normality
#> 1        5.48e-02        no
#> 2        1.58e-01        no
#> 3        1.10e-01        no
#> 4        8.66e-07        no
#> 5        2.73e-02        no
#> 6        8.70e-02        no
#> 7        4.60e-01        no
#> 8        4.65e-01        no
#> 9        2.58e-01        no
#> 10       2.72e-01        no
#> 11       3.38e-01        no
#> 12       1.81e-01        no
```
