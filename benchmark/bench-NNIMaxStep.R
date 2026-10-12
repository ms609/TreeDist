source("benchmark/_init.R")

Benchmark(NNIMaxStep(300))
Benchmark(NNIMaxStep(300, normalize = TRUE))
Benchmark(NNIMaxStep(2000, normalize = TRUE))
