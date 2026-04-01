# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla < memcheck/thisfile.R
devtools::build_vignettes(install = FALSE)
