## Test environments
* Local Windows 10 installation, R 4.0.2
* Windows 10, with `check_win_devel()`, R devel
* Ubuntu 16.04.6 LTS, R 3.6.0, release and devel, via
  [Travis CI](https://travis-ci.org/ms609/TreeDist)
* Mac OS X 10.13.6, R release, via Travis
* R-hub, with `check_for_cran()` and `check_with_sanitizers()`.

`check_with_sanitizers()` fails due to an error in the required package 
'phangorn', and I'm not aware of another way to verify that all ASAN issues
are fixed (as a Windows user).


## R CMD check results
There were no ERRORs or WARNINGs.

There were two NOTEs:
> Maintainer: 'Martin R. Smith <martin.smith@durham.ac.uk>'

> Suggests or Enhances not in mainstream repositories:
>   TreeDistData
> Availability using Additional_repositories specification:
>   TreeDistData   yes   https://ms609.github.io/packages/TreeDistData
[...]
> Package suggested but not available for checking: 'TreeDistData'

'TreeDistData' is too large to be submitted on CRAN.

All calls to `data(package = 'TreeDistData')` are wrapped within 
`if(require('TreeDistData')` to ensure that vignettes fail gracefully when
'TreeDistData' is not installed.


## Downstream dependencies
There is one downstream dependency, 'TBRDist' (which I maintain).

No changes to worse were identified by a local R CMD CHECK.
