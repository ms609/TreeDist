release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you updated REFERENCES.bib with a citation to the published study?",
    "Have you updated inst/CITATION with a citation to the published study?",
    "Have you updated the version number in inst/CITATION, NEWS and DESCRIPTION?",
    "Have you updated .zenodo.json?"
  )
}

# Additional tests:
# 
# check_win_devel(); check_rhub()
# revdepcheck()
# build_vignettes()
# 