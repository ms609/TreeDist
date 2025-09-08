library("TreeTools")

test_that("is.HPart() succeeds", {
  expect_true(is.HPart(as.HPart(TreeTools::BalancedTree(7))))
  expect_true(is.HPart(structure(class = "HPart",
                                 list(list("t1"), list("t2", "t3")))))
  expect_false(is.HPart(structure(class = "HPart", 
                                  list("t1", list("t2", "t3")))))
  expect_false(is.HPart(structure(class = "NonPart", 
                                  list(list("t1"), list("t2", "t3")))))
})

test_that("ReplicateHPart()", {
  h <- as.HPart(BalancedTree(6))
  expect_equal(ReplicateHPart(h, setNames(paste0("T", 1:6), paste0("t", 1:6))),
               rapply(h, toupper, how = "replace"))
})

test_that("HMI results match hmi.pynb", {
  # Non-hierarchical
  p1 <- list(list(19, 18, 5), list(14, 16, 3), list(7), list(10, 8), list(1, 17, 9, 4, 6, 15), list(2, 13, 11), list(12, 0))
  p2 <- list( list(12, 9), list(4, 2, 0, 7), list(16), list(5), list(8, 3, 1, 14), list(11, 6, 10), list(18, 17, 19), list(13, 15))
  expect_equal(HMIR(p1, p2), c(20, 0.9410980357245466))
  
  # Hierarchical
  hp1 <- list(list(23), list(list(list(list(list(list(16),  list(17)))))), list(list(12),  list(22, 13)), list(5), list(7), list(24), list(list(list(9),  list(list(14, 2))),  list(list(list(list(list(list(27),  list(3))))))), list(20, 29, 18), list(4), list(26, 15), list(list(10),  list(21, 25)), list(11), list(list(0, 28),  list(1),  list(6)), list(19, 8))
  hp2 <- list(list(list(list(0, 25),  list(24)),  list(6),  list(11, 28),  list(8)), list(list(list(19),  list(list(list(list(21),  list(4),  list(list(list(list(list(22, 7))))))))),  list(5)), list(list(3),  list(10, 23, 14)), list(list(27, 1, 16, 13, 18, 26, 9),  list(list(list(list(15),  list(list(list(list(list(list(12, 17)))))))),  list(2, 20)),  list(29)))
  
  expect_equal(HMIR(hp1, hp2), c(30, 1.0591260408329395))
})

test_that("HMI calculated correctly", {
  bal6 <- BalancedTree(6)
  pec6 <- PectinateTree(6)
  expect_equal(HierarchicalMutualInfo(bal6, pec6),
               HierachicalMutual(bal6, pec6))
})
