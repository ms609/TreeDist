## Test environments
* local Windows 10 install, R 3.4.1
* ubuntu 12.04 (on travis-ci), R 3.4.0 and devel

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Smith Martin R. <martins@gmail.com>'
  
  New submission
  
This package was previously submitted, and reviewed by Swetlana Herbrandt, who identified 
issues that we have now resolved:

- Copyright holders have now been added to the Authors@R field
- A reference has been added to the Descroption field of the DESCRIPTION file.

The previous submission contained a further note, which I understood to be a false positive.
This note is no longer showing on my checks:
* checking compiled code ... NOTE
  File 'TreeSearch/libs/x64/TreeSearch.dll':
  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'

  It is good practice to register native routines and to disable symbol
  search.
  
  See 'Writing portable packages' in the 'Writing R Extensions' manual.

## Downstream dependencies
There are currently no downstream dependencies for this package.