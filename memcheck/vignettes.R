# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla < memcheck/thisfile.R
pkgdown::build_articles(preview = FALSE)
