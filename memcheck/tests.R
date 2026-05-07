# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla < memcheck/thisfile.R
testthat::test_local()
