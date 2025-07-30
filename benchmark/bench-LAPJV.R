ub <- microbenchmark::microbenchmark
set.seed(1)
UnifMat <- function(n) matrix(runif(n * n), n, n)
test40 <- UnifMat(40)
test400 <- UnifMat(400)
test2000 <- UnifMat(2000)

ub(LAPJV(test40), LAPJV(test400), times = 100)
ub(LAPJV(test2000), times = 25)
