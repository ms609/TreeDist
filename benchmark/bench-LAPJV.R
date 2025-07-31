ub <- microbenchmark::microbenchmark
library("TreeDist") # not load_all
set.seed(1)
UnifMat <- function(n) matrix(runif(n * n), n, n)
test40 <- UnifMat(40)
test400 <- UnifMat(400)
test2000 <- UnifMat(1999) # non-2 multiple

ub(LAPJV(test40), LAPJV(test400), times = 100)
# With stl vectors in caller only
# 

ub(LAPJV(test2000), times = 25)
# With stl vectors in caller only
# Unit: milliseconds

