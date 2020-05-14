context("median.R")
library('TreeTools')
test_that("Median is calculated", {
  tenTrees <- as.phylo(1:10, 8L)
  expect_equal(as.phylo(4, 8), median(tenTrees))
  expect_equal(as.phylo(3, 8),
               median(tenTrees, Distance = RobinsonFoulds, breakTies = TRUE))
  expect_equivalent(as.phylo(c(3:5, 7), 8),
                    median(tenTrees, Distance = RobinsonFoulds, 
                           breakTies = FALSE))
})