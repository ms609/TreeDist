## Test environments
* local Windows 10 install, R 3.4.1
* ubuntu 12.04 (on travis-ci), R 3.4.2

## R CMD check results
 
I believe the NOTE may be a false positive - if not, I can't work out how to get around it.
I've followed the instructions at 
https://cran.r-project.org/doc/manuals/R-exts.html#Registering-native-routines
and the suggestions at
https://stackoverflow.com/questions/42313373/r-cmd-check-note-found-no-calls-to-r-registerroutines-r-usedynamicsymbols
to no avail.

## Downstream dependencies
There are currently no downstream dependencies for this package.