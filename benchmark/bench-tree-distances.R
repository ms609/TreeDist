source("benchmark/_init.R")

tr50 <- as.phylo(0:99, 50)
tr180 <- as.phylo(0:39, 180)
tr280 <- as.phylo(0:49, 280)

Benchmark(ClusteringInfoDistance(tr50))
Benchmark(ClusteringInfoDistance(tr280))
Benchmark(PhylogeneticInfoDistance(tr50))
Benchmark(PhylogeneticInfoDistance(tr180))

Benchmark(RobinsonFoulds(tr50))
Benchmark(RobinsonFoulds(tr280))
Benchmark(RobinsonFoulds(tr280), min_time = 20)

Benchmark(MutualClusteringInfo(tr50))
Benchmark(MutualClusteringInfo(tr280))
