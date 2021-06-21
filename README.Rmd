---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# `ccdf`

[![Downloads](https://cranlogs.r-pkg.org/badges/ccdf?color=blue)](https://www.r-pkg.org/pkg/ccdf)

## Overview

`ccdf` is a package for performing single-cell RNA-seq differential expression analysis and more generally **complex hypothesis testing**. 

The main function of the package is `ccdf_testing()`. It  allows to use either an asymptotic test for large sample size or a permutation test for small sample size with the argument `method`.

The methods implemented in this package are detailed in the following article:

> Gauthier M, Agniel D, Thiébaut R & Hejblum BP (2020). Distribution-free complex hypothesis testing for single-cell RNA-seq differential expression analysis, *BioRxiv*  [DOI:10.1101/2021.05.21.445165]( https://doi.org/10.1101/2021.05.21.445165)



## Installation

**To install `ccdf`, you can download the development version on [GitHub](https://github.com/Mgauth/ccdf).**
```{r GitHub-install, eval=FALSE}
#install.packages("devtoos")
devtools::install_github("Mgauth/ccdf")
```


## Example

Here is a basic example which shows how to use `ccdf` with simple generated data.

```{r example, eval = FALSE}
## Data Generation
X <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
Y <- replicate(10, ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,0.5,1)))
Y <- t(Y)
```

```{r estimation, eval = FALSE}
# Hypothesis testing
res_asymp <- ccdf_testing(exprmat=Y, variable2test=X, test="asymptotic") # asymptotic test
res_perm <- ccdf_testing(exprmat=Y, variable2test=X, test="permutations",
                         adaptive=TRUE) # adaptive permutation test
```


 -- Marine Gauthier, Denis Agniel, Rodolphe Thiébaut & Boris Hejblum
