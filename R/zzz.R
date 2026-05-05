.onUnload <- function(libpath) {
  StopParallel(quietly = TRUE)
  library.dynam.unload("TreeDist", libpath)
}

## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you updated inst/REFERENCES.bib with a full citation to Smith & Donoghue 2022?"
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
# # tools::resaveRdaFiles("R", compress="auto") - is default bzip2 the optimal?
# # tools::checkRdaFiles("R") - set optimal compression in `data-raw`
