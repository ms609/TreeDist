context("median.R")
library('TreeTools')
test_that("Median is calculated", {
  tenTrees <- as.phylo(1:10, 8L)
  expect_equal(4, median(tenTrees, index = TRUE))
  expect_equal(as.phylo(4, 8), median(tenTrees))
  expect_equal(as.phylo(3, 8),
               median(tenTrees, Distance = RobinsonFoulds, breakTies = TRUE))
  expect_equal(c(3:5, 7), median(tenTrees, Distance = RobinsonFoulds, 
                                 index = TRUE, breakTies = FALSE))
  expect_equivalent(as.phylo(c(3:5, 7), 8),
                    median(tenTrees, Distance = RobinsonFoulds, 
                           breakTies = FALSE))
})