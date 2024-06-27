library("TreeDist")
t1 <- Postorder(as.phylo(0:5, 6))
t2 <- Postorder(BalancedTree(6))
t3 <- Postorder(PectinateTree(6))
expect_equal(PathDist(UnrootTree(t1), UnrootTree(t2)),
             phangorn::path.dist(t1, t2))
expect_equal(PathDist(UnrootTree(t2), UnrootTree(t3)),
             phangorn::path.dist(t3, t2))

postTrees <- Postorder(as.phylo(0:5, 182))
expect_equal(PathDist(postTrees), phangorn::path.dist(postTrees))


ub <- microbenchmark::microbenchmark
ub(PathDist(postTrees), phangorn::path.dist(postTrees))
pv <- profvis::profvis
pv(ub(PathDist(postTrees), phangorn::path.dist(postTrees)))
