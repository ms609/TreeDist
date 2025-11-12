# CONTRIBUTING

### Fixing typos

Small typos or grammatical errors in documentation may be edited
directly using the GitHub web interface, so long as the changes are made
in the *source* file.

- DO: edit a roxygen comment in a `.R` file below `R/`.
- DON’T: edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an
issue and make sure someone from the team agrees that it’s a problem. If
you’ve found a bug, create an associated issue and illustrate the bug
with a minimal [reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

- We recommend that you create a Git branch for each pull request
  (PR).  
- Look at the GitHub actions and CodeCovr build status before and after
  making changes.
- We follow [Google’s R style
  guide](https://google.github.io/styleguide/Rguide.html)
- We use camelCase for variable names, and TitleCase for function names.
- We use the Oxford ending of ‘ize’ (not ‘ise’), and UK spelling
  (e.g. ‘colour’) where it is not possible to avoid the distinction
  (e.g. by shortening to ‘col’)
- We use [roxygen2](https://cran.r-project.org/package=roxygen2) to
  build documentation.
- We use [testthat](https://cran.r-project.org/package=testthat).
  Contributions with test cases included are easier to accept.
- We use [semantic versioning](https://semver.org/).
- For user-facing changes, add a bullet to the top of `NEWS.md` below
  the current development version header describing the changes made
  followed by your GitHub username, and links to relevant
  issue(s)/PR(s).

### Contributor license agreement

Contributors agree to reassign the copyright of their contributions to
the maintainers of the package.

### Code of Conduct

Please note that the project is released with a [Contributor Code of
Conduct](https://ms609.github.io/TreeDist/CODE_OF_CONDUCT.md). By
contributing to this project you agree to abide by its terms.

### Prefer to e-mail?

E-mail the person listed as maintainer in the `DESCRIPTION` file of this
repo. Private discussions over email don’t help others - but of course
email is totally warranted for any sensitive problems.

### Thanks for contributing!

This document is adapted from the [tidyverse contributing
guide](https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md).
