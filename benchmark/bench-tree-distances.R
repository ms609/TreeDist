source("benchmark/_init.R")
library("TreeTools", quietly = TRUE)
library("TreeDist")

tr50 <- as.phylo(0:99, 50)
tr200 <- as.phylo(0:39, 200)

Benchmark("CID50", ub(ClusteringInfoDistance(tr50), times = 20))
Benchmark("CID200", ub(ClusteringInfoDistance(tr200), times = 10))
Benchmark("PID50", ub(PhylogeneticInfoDistance(tr50), times = 20))
Benchmark("PID200", ub(PhylogeneticInfoDistance(tr200), times = 10))
Benchmark("RF50", ub(RobinsonFoulds(tr50)))
Benchmark("RF200", ub(RobinsonFoulds(tr200)))

Benchmark("MuClInf50", ub(MutualClusteringInfo(tr50), times = 20))
Benchmark("MuClInf200", ub(MutualClusteringInfo(tr200), times = 20))
