source("benchmark/_init.R")

UnifMat <- function(n) matrix(runif(n * n), n, n)
test40 <- UnifMat(40)
test400 <- UnifMat(400)
test2000 <- UnifMat(1999) # non-2 multiple

Benchmark("LAPJV40", ub(LAPJV(test40), LAPJV(test400), times = 100))
Benchmark("LAPJV40", ub(LAPJV(test2000), times = 25))
