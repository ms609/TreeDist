source("benchmark/_init.R")

UnifMat <- function(n) matrix(runif(n * n), n, n)
test40 <- UnifMat(40)
test400 <- UnifMat(400)
test2000 <- UnifMat(1999) # non-2 multiple

Benchmark(LAPJV(test40))
Benchmark(LAPJV(test400))
Benchmark(LAPJV(test2000))
