library('TreeTools')
context('tree_distance_utilities.R [Day]')

test_that("DayDistance() is robust", {
  #TODO lots!
  
  PrepareTree <- function (tree) {
    RenumberTips(tree, as.character(seq_along(tree$tip.label)))
  }
  PrepareText <- function (text) {
    PrepareTree(ape::read.tree(text = text))
  }
  
  TestRF <- function (t1, t2) {
    expect_equal(RobinsonFoulds(t1, t2), 
                 DayDistance(DayRobinsonFoulds,
                             as.ClusterTable(t1),
                             as.ClusterTable(t2))[1])
  }
  
  t1 <- PrepareText("(1, (2, (3, (4, (5, 6)))));")
  t2 <- PrepareText("(1, (2, ((4, 3), (6, 5))));")
  t3 <- PrepareText("(1, (((2, 3), 5), (4, 6)));")
  
  TestRF(t1, t1)
  TestRF(t1, t2)
  TestRF(t2, t1)
  TestRF(t1, t3)
  TestRF(t3, t2)
  
  clist <- list(as.ClusterTable(t1),
                as.ClusterTable(t2),
                as.ClusterTable(t3))
  
  expect_equivalent(RobinsonFoulds(list(t1, t2, t3)),
                    6L - robinson_foulds_all_pairs(clist))
  expect_equivalent(RobinsonFoulds(list(t1, t2, t3)),
                    DayDistance(DayRobinsonFoulds, clist, clist))
  
  t1p <- ape::read.tree(text = "(1, (2, (3, (4, (5, 6)))));")
  t2p <- ape::read.tree(text = "(1, (2, ((4, 3), (6, 5))));")
  t3p <- ape::read.tree(text = "(1, (((2, 3), 5), (4, 6)));")
  microbenchmark::microbenchmark(phangorn::RF.dist(c(t1p, t2p, t3p)),
                                 RobinsonFoulds(list(t1, t2, t3)),
                                 DayDistance(DayRobinsonFoulds, list(as.ClusterTable(t1),
                                                                     as.ClusterTable(t2),
                                                                     as.ClusterTable(t3)), clist),
                                 DayDistance(DayRobinsonFoulds, clist, clist))

  t1p <- BalancedTree(seq_len(100))
  t2p <- PectinateTree(seq_len(100))
  t3p <- LeafLabelInterchange(t1p, 3)
  t4p <- LeafLabelInterchange(t2p, 3)
  t5p <- LeafLabelInterchange(t1p, 8)
  pTr <- c(t1p, t2p, t3p, t4p, t5p)
  prTr <- lapply(pTr, PrepareTree)
  cls <- lapply(prTr, as.ClusterTable)
  microbenchmark::microbenchmark(phangorn::RF.dist(pTr),
                                 RobinsonFoulds(pTr),
                                 DayDistance(DayRobinsonFoulds, lapply(prTr, as.ClusterTable),
                                                                cls),
                                 lapply(prTr, as.ClusterTable),
                                 robinson_foulds_all_pairs(cls),
                                 DayRobinsonFoulds[[1]](cls),
                                 DayDistance(DayRobinsonFoulds, cls, cls))
  cls <- c(cls, cls)
  cls <- c(cls, cls)
  cls <- c(cls, cls)
  cls <- c(cls, cls)
  microbenchmark::microbenchmark(DayDistance(DayRobinsonFoulds, cls, cls))
  profvis::profvis(DayDistance(DayRobinsonFoulds, cls, cls))
  
  
})
  