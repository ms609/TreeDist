[![Build Status](https://travis-ci.org/ms609/TreeSearch.svg?branch=master)](https://travis-ci.org/ms609/TreeSearch)
[![codecov](https://codecov.io/gh/ms609/TreeSearch/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/TreeSearch)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/TreeSearch)](https://cran.r-project.org/package=TreeSearch)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/TreeSearch)](https://cran.r-project.org/package=TreeSearch)<!--[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)
-->
[![Project Status: Active – – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

# TreeSearch
This package exists to allow parsimony-style tree searches in R.

It extends the basic functionality available in phangorn, with a view to making tree search faster and more efficient, 
and allowing user-defined optimality criteria to be employed.

Implied weighting and heuristic searches such as the Parsimony Ratchet are implemented (function: `Ratchet`).
# Installation

Install and load the library from CRAN as follows:
```
install.packages('TreeSearch')
library('TreeSearch')
```

If you're feeling brave, you can install the development version thus:
```r
if(!require(devtools)) install.packages("devtools")
devtools::install_github('ms609/TreeSearch')
```


# Optimality criteria

TreeSearch allows the implementation of various optimality criteria, including
- The Profile Parsimony approach introduced by Faith and Trueman (2001).
- Successive Approximations weighting (Farris 1969).
- Fitch parsimony with inapplicable data (Brazeau, Guillerme and Smith, 2019).

It is also possible to specify bespoke optimality criteria.


# Documentation

- [Package functions](https://cran.r-project.org/web/packages/TreeSearch/TreeSearch.pdf) reference manual
- [Getting started](https://cran.r-project.org/web/packages/TreeSearch/vignettes/getting-started.html)
- [Loading phylogenetic data into R](https://cran.r-project.org/web/packages/TreeSearch/vignettes/load-data.html)
- [Parsimony search with inapplicable data](https://cran.r-project.org/web/packages/TreeSearch/vignettes/inapplicable.html)

- [Calculating concavity profiles](https://cran.r-project.org/web/packages/TreeSearch/vignettes/profile-scores.html) for Profile Parsimony
- [Tree search with profile parsimony](https://cran.r-project.org/web/packages/TreeSearch/vignettes/profile.html)


# References

Brazeau, M. D., Guillerme, T. and Smith, M. R. in press. 
  An algorithm for morphological phylogenetic analysis with inapplicable data. 
  Systematic Biology. doi:[10.1093/sysbio/syy083](https://dx.doi.org/10.1093/sysbio/syy083).

D. P. Faith, J. W. H. Trueman, Towards an inclusive philosophy for phylogenetic inference.
  Syst. Biol. 50, 331–350 (2001).   doi:[10.1080/10635150118627](https://dx.doi.org/10.1080/10635150118627)

Farris, J. S. (1969). A successive approximations approach to character weighting. 
  Systematic Biology, 18(4), 374–385.  doi:[10.2307/2412182](https://dx.doi.org/10.2307/2412182)
