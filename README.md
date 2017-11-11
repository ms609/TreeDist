[![Build Status](https://travis-ci.org/ms609/TreeSearch.svg?branch=master)](https://travis-ci.org/ms609/TreeSearch)
[![codecov](https://codecov.io/gh/ms609/TreeSearch/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/TreeSearch)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/TreeSearch)](https://cran.r-project.org/package=TreeSearch)
<!--[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/TreeSearch)](https://cran.r-project.org/package=TreeSearch)-->
<!--[![Research software impact](http://depsy.org/api/package/cran/TreeSearch/badge.svg)](http://depsy.org/package/r/TreeSearch)-->
[![Project Status: Active – – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

# TreeSearch
This package exists to allow parsiomony-style tree searches in R.

It extends the basic functionality available in phangorn, with a view to making tree search faster and more efficient, 
and allowing user-defined optimality criteria to be employed.

Heuristic searches such as the Parsimony Ratchet are implemented (function: `Ratchet`).

The key function is `TreeSearch`, which takes a tree and a dataset; functions can be specified to 'load' the 
data (perhaps sending it to C?) and to score a tree (the Fitch algorithm is used by default).

# Installation

Install and load the library from CRAN as follows:
```
install.packages('TreeSearch')
library('TreeSearch')
```

The library requires a working version of phangorn > 2.2.1.  The version on the CRAN repository at 2 Nov 2017 
has caused some issues during installation; an alternative is to install from a known working version of 30 Oct 2017:

```
if (!require(devtools)) install.packages('devtools')
devtools::install_github('KlausVigo/phangorn', ref='7192bfb4403c35c16a7b735160525d272736b061') 
```
