.onUnload <- function (libpath) {
  library.dynam.unload("TreeDist", libpath)
}

## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you checked the Vignettes for sanity?",
    "Have you checked pkgdown::build_reference_index()?",
    "Have you updated README.md, inst/REFERENCES.bib & inst/CITATION with a citation to the published study?",
    "Have you updated the version number in .zenodo.json, NEWS & DESCRIPTION?"
    )
}

# Additional tests:
# 
# build_reference_index()
# codemetar::write_codemeta()
# check_win_devel(); rhub::check_for_cran()
# revdepcheck::revdep_check()
# build_vignettes()
# tools::resaveRdaFiles('R', compress='auto') - is default bzip2 the optimal?
# tools::checkRdaFiles('R') - set optimal compression in `data-raw`
