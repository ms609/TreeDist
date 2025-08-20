source("benchmark/_init.R")

tr50 <- as.phylo(0:99, 50)
tr200 <- as.phylo(0:39, 200)

Benchmark(ClusteringInfoDistance(tr50))
Benchmark(ClusteringInfoDistance(tr200))
Benchmark(PhylogeneticInfoDistance(tr50))
Benchmark(PhylogeneticInfoDistance(tr200))

Benchmark(RobinsonFoulds(tr50))
Benchmark(RobinsonFoulds(tr200))

Benchmark(MutualClusteringInfo(tr50))
Benchmark(MutualClusteringInfo(tr200))
