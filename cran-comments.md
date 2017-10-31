## Test environments
* local Windows 10 install, R 3.4.1
* ubuntu 12.04 (on travis-ci), R 3.4.0 and devel

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Smith Martin R. <martins@gmail.com>'
  
  New submission
  
This is indeed a new submission.

* checking compiled code ... NOTE
  File 'TreeSearch/libs/x64/TreeSearch.dll':
  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'

  It is good practice to register native routines and to disable symbol
  search.
  
  See 'Writing portable packages' in the 'Writing R Extensions' manual.

I've followed the instructions at 
https://cran.r-project.org/doc/manuals/R-exts.html#Registering-native-routines
and the suggestions at
https://stackoverflow.com/questions/42313373
to no avail.  As such, I suspect that the latter NOTE may be a false positive - if not, I can't work out how to get around it.

## Downstream dependencies
There are currently no downstream dependencies for this package.