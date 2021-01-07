## Test environments
* Microsoft Windows:
  * Local Windows 10 installation, R 4.0.3
  * windows-latest: Microsoft Windows Server 2019, Windows 10.0.17763, 
    R release, via Github Actions
  * win_devel: with `check_win_devel()`, R devel
  
* Linux:
  * Ubuntu 20.04.1 LTS, R 3.6.0, release and devel, via Github Actions
  * via R-hub, with `rhub::check_for_cran()` and `rhub::check_with_sanitizers()`.
  
* Mac OS X 10.15.7, R release, via Github actions


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
