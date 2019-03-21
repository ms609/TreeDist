## Test environments
### Windows 10:
* local Windows 10 install, R 3.5.2
* Windows 10 via check_win_devel(), R devel

### Linux:
* ubuntu 14.04.5 LTS (on travis-ci), R 3.4.0 and release
* Using check_rhub()

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

Kurt Hornik requested that I resubmit this package immediately in order
to prepare for the new random number generation system in R 3.6.0.

## Downstream dependencies
There are currently two downstream dependencies for this package.

`revdepcheck::revdep_check()` does not yet support the Windows platform.
`devtools::revdep_check()` reported no packages with problems.
