# TreeSearch
This package exists to allow parsiomony-style tree searches in R.

It extends the basic functionality available in phangorn, with a view to making tree search faster and more efficient.
Heuristic searches such as the Parsimony Ratchet are implemented.

The library is currently in development and is not ready for use.

The library requires the latest versions of ape and phangorn:

```
install.packages('ape', '.', repos = 'http://ape-package.ird.fr/')
devtools::install_github('KlausVigo/phangorn', ref='1167f0be62f13cfad0fca8ae8224318c407195bf')
```