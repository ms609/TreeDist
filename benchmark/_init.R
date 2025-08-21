library("TreeTools")
set.seed(54321)

Benchmark <- function(..., min_iterations = NULL) {
  result <- bench::mark(..., min_iterations = min_iterations %||% 3,
                        time_unit = "us")
  if (interactive()) {
    print(result)
  } else {
    fileroot <- gsub("[\"']", "",
                     gsub("[\\(\\):, /]", "_", as.character(result$expression)))
    .FileName <- function(fileRoot, i) {
      paste0(c(fileroot, i, "bench.Rds"), collapse = ".")
    }
    i <- double(0)
    while(file.exists(.FileName(fileroot, i))) {
      if (length(i) == 0) i <- 0
      i <- 1 + i
    }
    saveRDS(result, .FileName(fileroot, i))
  }
}
