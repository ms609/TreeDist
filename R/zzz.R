.onUnload <- function (libpath) {
  library.dynam.unload("TreeDist", libpath)
}

## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you checked the Vignettes for sanity?",
    "Have you checked pkgdown::build_reference_index()?",
    "Have you refreshed the package meta with codemetar::write_codemeta()?",
    "Have you updated README.md, inst/REFERENCES.bib & inst/CITATION with a citation to the published study?",
    "Have you updated the version number in inst/CITATION, .zenodo.json, NEWS & DESCRIPTION?"
    )
}

# Additional tests:
# 
# check_win_devel(); check_rhub()
# revdepcheck::revdep_check()
# build_vignettes()
# tools::resaveRdaFiles('data', compress='auto' - is default bzip2 the optimal?
# tools::checkRdaFiles('data') - set optimal compression in `data-raw`
