source("benchmark/_init.R")

Benchmark(NNIMaxStep(300))
Benchmark(NNIMaxStep(300, normalize = TRUE))
