## Test environments
* Local Windows 10 installation, R 3.6.1
* Windows 10, with `check_win_devel(quiet = TRUE)`, R devel
* Ubuntu 16.04.6 LTS, R 3.5.0, release and devel, via [Travis CI](https://travis-ci.org/ms609/TreeDist)
* Mac OS X 10.13.6, R release, via Travis
* R-hub, with `check_for_cran()`

## R CMD check results
There were no ERRORs or WARNINGs.

There were two NOTEs:
> Maintainer: 'Martin R. Smith <martin.smith@durham.ac.uk>'
> 
> New submission

This is a new submission.

> Possibly mis-spelled words in DESCRIPTION:
>    Bogdanowicz (19:32)
>    Bocker (17:38)
>    Colijn (21:15)
>    Foulds (15:42, 17:20)
>    Giaro (19:46)
>    Jaccard (17:3)
>    NNI (22:34)
>    Nye (18:17)
>    al (17:48, 18:24, 22:72)
>    et (17:45, 18:21, 22:69)

These are the content of references, except the defined acronym NNI (which is
more familiar to some users than its spelled-out version).

> * checking for unstated dependencies in vignettes ... NOTE
> 'library' or 'require' call not declared from: 'TreeDistData'
> 'data(package=)' call not declared from: 'TreeDistData'

'TreeDistData', which depends on 'TreeDist', is not yet available on CRAN, 
so cannot be added to the Suggests field in the DESCRIPTION file.

All calls to `data(package='TreeDistData')` are wrapped within 
`if(require('TreeDistData')` to ensure that vignettes fail gracefully if
'TreeDistData' is not installed.


## Downstream dependencies
There are currently no downstream dependencies for this package.
