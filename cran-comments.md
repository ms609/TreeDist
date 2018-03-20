## Test environments
### Windows 10:
* local Windows 10 install, R 3.4.3
* Windows x86_64-w64-mingw32, via R-hub, R 3.4.4
* Windows 10 via win_build(), R devel

### Linux:
* ubuntu 14.04.5 (on travis-ci), R 3.4.0 and release
* Debian, Ubuntu, Fedora and Centos, via R-hub, R 3.4.3

## R CMD check results
There were no ERRORs or WARNINGs. 

There was one NOTE:

> NOTE
> Maintainer: 'Martin R. Smith <martins@gmail.com>'
> 
> Days since last update: 5

Version 0.1.0 fixed an error in the C code that had been flagged by the CRAN package check.  Fixing this error exposed some new C warnings, which this version 0.1.1 addressed.

Uwe pointed out an issue building vignettes on certain platforms, and I've removed the probelmatic code in the present submission, v0.1.2.


## Downstream dependencies
There are currently no downstream dependencies for this package.
