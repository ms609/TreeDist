# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla < memcheck/thisfile.R
library("TreeDist")

example_db <- tools::Rd_db("TreeDist")

cat("Running examples for", length(example_db), "topics\n")

failures <- character(0)

for (topic in names(example_db)) {
  cat("\n>>> Example:", topic, "\n")
  ex <- tools::Rd2ex(example_db[[topic]])
  
  if (length(ex) == 0L) {
    cat("No example found for topic:", topic, "\n")
    next
  }
  
  ex_file <- tempfile(fileext = ".R")
  writeLines(ex, ex_file)
  
  # Try running example code in globalenv, catching errors
  tryCatch(
    {
      sys.source(ex_file, envir = globalenv())
      cat("\U2713 Success:", topic, "\n")
    },
    error = function(e) {
      cat("\U2718 Error in topic:", topic, "\n", conditionMessage(e), "\n")
      failures <<- c(failures, topic)
    }
  )
}
cat("\nFinished running examples.\n")

if (length(failures)) {
  cat("\U274c Failures in", length(failures), "topics:\n")
  print(failures)
  quit(status = 1)
} else {
  cat("\U2705 All examples ran successfully.\n")
}
