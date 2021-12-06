## Test environments
* Microsoft Windows:
  * Local Windows 10 installation, R 4.1.2
  * windows-latest: Microsoft Windows Server 2019, Windows 10.0.17763, 
    R release, via [GitHub Actions](https://github.com/ms609/TreeDist/actions)
  * win_devel: with `devtools::check_win_devel()`, R devel
  * win_oldrel: with `devtools::check_win_oldrelease()`.
  
* Linux:
  * Ubuntu 20.04.1 LTS, R 4.1, release and devel, via GitHub Actions
  * [valgrind mem-check](https://github.com/ms609/TreeDist/actions/workflows/memcheck.yml)
  * via R-hub, with `rhub::check_for_cran()`.
  
* Mac OS X 10.15.7, R release, via GitHub actions


## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE:

> Suggests or Enhances not in mainstream repositories:
>   TreeDistData
> Availability using Additional_repositories specification:
>   TreeDistData   yes   https://ms609.github.io/packages/TreeDistData
[...]
> Package suggested but not available for checking: 'TreeDistData'

'TreeDistData' is too large to be submitted on CRAN.

All calls to `data(package = 'TreeDistData')` are wrapped within 
`if(require('TreeDistData'))` to ensure that vignettes fail gracefully when
'TreeDistData' is not installed.


## Downstream dependencies
There is one downstream dependency, 'TBRDist' (which I maintain).

No changes to worse were identified by
["revdepcheck"](https://github.com/ms609/TreeDist/actions/workflows/revdepcheck.yml).
