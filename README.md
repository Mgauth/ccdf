
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `ccdf`

[![Downloads](https://cranlogs.r-pkg.org/badges/ccdf?color=blue)](https://www.r-pkg.org/pkg/ccdf)

## Overview

`ccdf` is a package for performing single-cell RNA-seq differential
expression analysis and more generally ***complex hypothesis testing***.

The main function of the package is `ccdf_testing()`. It allows to use
either an asymptotic test for large sample size or a permutation test
for small sample size with the argument
`method`.

<!-- The method implemented in this package is detailed in the following article: -->

<!-- > Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2020). ..., *bioRxiv* ... . [DOI: .../...](url) -->

## Installation

***To install `ccdf`, you can download the development version on
[GitHub](https://github.com/Mgauth/ccdf)***

``` r
#install.packages("devtoos")
devtools::install_github("Mgauth/ccdf")
```

– Marine Gauthier, Denis Agniel, Rodolphe Thiébaut & Boris Hejblum
