[![Build Status](https://travis-ci.org/ms609/TreeSearch.svg?branch=master)](https://travis-ci.org/ms609/TreeSearch)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

# TreeSearch
This package exists to allow parsiomony-style tree searches in R.

It extends the basic functionality available in phangorn, with a view to making tree search faster and more efficient.
Heuristic searches such as the Parsimony Ratchet are implemented.

The library requires the latest (working) versions of ape and phangorn:

```
install.packages('ape', '.', repos = 'http://ape-package.ird.fr/')
devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
```
