[![Build Status](https://travis-ci.org/ms609/TreeSearch.svg?branch=master)](https://travis-ci.org/ms609/TreeSearch)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

# TreeSearch
This package exists to allow parsiomony-style tree searches in R.

It extends the basic functionality available in phangorn, with a view to making tree search faster and more efficient, 
and allowing user-defined optimality criteria to be employed.

Heuristic searches such as the Parsimony Ratchet are implemented (function: `Ratchet`).

The library requires a working version of phangorn > 2.2.1, which is not yet available through CRAN:

```
devtools::install_github('KlausVigo/phangorn', ref='7192bfb4403c35c16a7b735160525d272736b061') # 30 oct 2017
```

The key function is `TreeSearch`, which takes a tree and a dataset; functions can be specified to 'load' the 
data (perhaps sending it to C?) and to score a tree (the Fitch algorithm is used by default).