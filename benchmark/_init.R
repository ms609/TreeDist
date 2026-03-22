library("TreeTools", quietly = TRUE)
library("TreeDist")
ub <- microbenchmark::microbenchmark
set.seed(1337)

Benchmark <- function(id, result) {
  if (interactive()) {
    print(result)
  } else {
    saveRDS(result, paste0(id, ".bench.Rds"))
  }
}
