## Test environments
* local Windows 10 install, R 3.4.3
* ubuntu 14.04.5 (on travis-ci), R 3.4.0 and release

## R CMD check results
There were no ERRORs or WARNINGs. 

There was one NOTE:

> NOTE
> Maintainer: 'Martin R. Smith <martins@gmail.com>'
> 
> Days since last update: 1

Version 0.1.0 fixed an error in the C code that had been flagged by the CRAN package check.  Fixing this error exposed some new warnings, which this resubmission addresses.

I'm sorry that this has meant multiple submissions in a short space of time.  If there is a way for me to test the package in the environments that CRAN uses ahead of submission, to avoid this situation recurring, please do let me know.


## Downstream dependencies
There are currently no downstream dependencies for this package.
