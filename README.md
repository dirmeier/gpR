<h1 align="center"> lvgpp </h1>

[![Build Status](https://travis-ci.org/rafstraumur/lvgpp.svg?branch=master)](https://travis-ci.org/rafstraumur/lvgpp.svg?branch=master)
[![codecov](https://codecov.io/gh/rafstraumur/lvgpp/branch/master/graph/badge.svg)](https://codecov.io/gh/rafstraumur/lvgpp)


Gaussian processes for machine learning in R.

## Introduction

Gaussian Processes have recently gained a lot of attention in machine learning. Since I am fan of prediction using GPs as well I'll introduce the theory roughly. The package <code>lvgpp</code> shows how training (calculation of the posterior predictive) and prediction is done when the kernel parameters are known. In the next versions I will also implement how those are calculated by optimizing marginal likelihood and probably include more kernels.

## Installation
 
Install `lvgpp` using:

```{r}
library(devtools)
install_github("rafstraumur/lvgpp") 
```

from the R-console.

## Usage

Load the package using `library(lvgpp)`. We provide a vignette for the package that can be called using: `vignette("lvgpp")`. This should be all the information you need. For regression try the demo-tour using:

```{r}
demo.regression()
```

or for classification (i.e. binomial responses):

```{r}
demo.bin.classification()
```

Also check out the source code for more info or just write me!

## Author

* Simon Dirmeier <a href="mailto:simon.dirmeier@gmx.de">simon.dirmeier@gmx.de</a>