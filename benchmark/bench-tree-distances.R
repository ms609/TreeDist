library("TreeTools", quietly = TRUE)
library("TreeDist")
ub <- microbenchmark::microbenchmark

tr50 <- as.phylo(0:99, 50)
tr200 <- as.phylo(0:39, 200)
ub(times = 10,
   ClusteringInfoDistance(tr200), ClusteringInfoDistance(tr50),
   MatchingSplitDistance(tr200), MatchingSplitDistance(tr50)
   )
ub(RobinsonFoulds(tr200), RobinsonFoulds(tr50))
