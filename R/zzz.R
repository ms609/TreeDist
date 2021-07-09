.onUnload <- function (libpath) {
  library.dynam.unload("TreeDist", libpath)
}

## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you updated README.md, inst/REFERENCES.bib, man/KendallColijn & inst/CITATION with a citation to the published studies?"
    )
}

# Additional steps:
#
# Propogate changes in README.md to R/TreeDist-package.R

# Additional tests:
# 
# spell_check()
# pkgdown::build_reference_index()
#
# run_examples()
# build_vignettes()
# build_manual() # PDF support for special characters
#
# devtools::check_win_devel(quiet = TRUE); rhub::check_for_cran()
# Check valgrind section of GitHub actions for memcheck errors
# 
# revdepcheck::revdep_check()
#
# codemetar::write_codemeta()
#
# # Unnecessary:
# # tools::resaveRdaFiles('R', compress='auto') - is default bzip2 the optimal?
# # tools::checkRdaFiles('R') - set optimal compression in `data-raw`
