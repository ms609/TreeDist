source("benchmark/_init.R")

tr <- as.phylo(seq_len(1000), 183)
sp <- as.Splits(tr)

Benchmark(ClusteringInfoDistance(tr[[1]], tr))
Benchmark(ClusteringInfoDistance(sp[[1]], sp))
Benchmark(MutualClusteringInfo(tr[1:100], tr[1:100 + 500]))
