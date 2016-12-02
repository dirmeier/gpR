<h1 align="center"> gpR </h1>

[![Build Status](https://travis-ci.org/dirmeier/gpR.svg?branch=master)](https://travis-ci.org/dirmeier/gpR)
[![codecov](https://codecov.io/gh/dirmeier/gpR/branch/master/graph/badge.svg)](https://codecov.io/gh/dirmeier/gpR)

Gaussian processes for machine learning in R and FORTRAN.

## Introduction

Gaussian Processes have recently gained a lot of attention in machine learning. <code>gpR</code> shows how training (calculation of the posterior predictive) and prediction is done when the kernel parameters are *known*. In the next versions I will implement how those are calculated by optimizing the marginal likelihood and probably include more kernels.

## Installation
 
Install `gpR` using:

```{r}
devtools::install_github("dirmeier/gpR") 
```

from the R-console.

## Usage

Load the package using `library(gpR)`. We provide a vignette for the package that can be called using: `vignette("gpR")`. This should be all the information you need. For regression try the demo-tour using:

```{r}
demo.regression()
```

or for classification (i.e. binomial responses):

```{r}
demo.bin.classification()
```

Also check out the source code for more info, fork the package, or just write me!

## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@gmx.de">simon.dirmeier@gmx.de</a>
