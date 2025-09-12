source("benchmark/_init.R")

set.seed(1)
tips <- 64
chars <- 128
char <- matrix(sample(1:3, tips * chars, replace = TRUE), tips, chars)
tree <- TreeTools::RandomTree(tips)

Benchmark(apply(char, 2, CharAMI, tree))
