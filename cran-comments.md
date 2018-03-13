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

Brian Ripley kindly pointed out that the recently-uploaded v0.0.8 triggered errors in certain build environments.  This resubmission fixes, I hope, the warnings and notes.

If there is a way for me to test the package in the environments that CRAN
uses ahead of submission, to avoid this situation recurring, I'd be glad
to hear of it.


## Downstream dependencies
There are currently no downstream dependencies for this package.
