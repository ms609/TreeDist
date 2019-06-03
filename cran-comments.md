## Test environments
* Windows 10 on local machine, R 3.6.0
* Windows 10 via check_win_devel(), R devel
* ubuntu 14.04.5 LTS (on travis-ci), R 3.4.0 and release
* Using check_rhub()

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

The previous release, v0.3.0, triggered a ` noLD ` error on building the vignettes.
I have revised the vignettes to reduce their compile time, and do not encounter
a `noLD` error when building the vignettes locally.

## Downstream dependencies
There are currently two downstream dependencies for this package.

Neither `revdepcheck::revdep_check()` nor 
`devtools::revdep_check()` were able to check downstream dependencies for 
errors.
