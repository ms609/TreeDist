library("TreeTools")
tr1 <- BalancedTree(256)
tr2 <- RandomTree(256)
tr3 <- BalancedTree(1234)
tr4 <- RandomTree(1234)
tr5 <- BalancedTree(25)
tr6 <- PectinateTree(25)

ub(RobinsonFoulds(tr1, tr2),
   RobinsonFoulds(tr3, tr4),
   RobinsonFoulds(tr5, tr6),
   times = 1000)
