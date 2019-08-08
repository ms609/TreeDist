release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you updated REFERENCES.bib with a citation to the published study?",
    "Have you updated inst/CITATION with a citation to the published study?",
    "Have you updated the version number in inst/CITATION, NEWS and DESCRIPTION?"
  )
}

# Additional tests:
# 
# check_win_devel(); check_rhub()
# revdepcheck::revdep_check()
# build_vignettes()
# tools::resaveRdaFiles('data', compress='auto') - is default of bzip2 the optimal?
# tools::checkRdaFiles('data') - set optimal compression in `data-raw`
