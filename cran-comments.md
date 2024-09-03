## Test environments
* Local machine, Windows 10, R 4.4.1 (2024-06-14 ucrt)

* `devtools::check_win_devel()`

* [GitHub Actions](https://github.com/ms609/TreeDist/actions):
  - windows-latest, R release
  - Ubuntu-latest, R 3.6, release and devel
  - Mac OS X 10.15.7, R release, via GitHub actions
  
* [valgrind mem-check](https://github.com/ms609/TreeDist/actions/workflows/memcheck.yml)

* [R-hub](https://github.com/ms609/TreeDist/actions/workflows/rhub.yaml):
  - linux,macos,macos-arm64,windows,atlas,clang-asan,valgrind


## Downstream dependencies

["revdepcheck"](https://github.com/ms609/TreeDist/actions/workflows/revdepcheck.yml)
confirmed no changes to worse in the three downstream dependencies:

  'Rogue', 'TBRDist', and 'TreeSearch'
  
(all of which I maintain).
  

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

All calls to `data(package = "TreeDistData")` are wrapped within 
`if(require('TreeDistData'))` to ensure that vignettes fail gracefully when
'TreeDistData' is not installed.

