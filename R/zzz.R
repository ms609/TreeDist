.onUnload <- function(libpath) {
  StopParallel()
  library.dynam.unload("TreeDist", libpath)
}

## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you updated README.md, DESCRIPTION, inst/REFERENCES.bib, man/KendallColijn & inst/CITATION with a citation to the published studies?"
    )
}

# Additional steps:
#
# codemeta::write_codemeta()
# 
# Propagate changes in README.md to R/TreeDist-package.R

# Additional tests:
#
# run_examples()
# build_vignettes()
# build_manual() # PDF support for special characters

# # Unnecessary:
# # tools::resaveRdaFiles('R', compress='auto') - is default bzip2 the optimal?
# # tools::checkRdaFiles('R') - set optimal compression in `data-raw`
