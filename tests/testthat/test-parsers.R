context("Parser functions: TNT")
test_that("Nexus file can be parsed", {
  filename <- 'test-parse-nexus.nexus'
  read <- ReadCharacters(filename)
  expect_equal(192, ncol(read))
  expect_equal(80, nrow(read))
  expect_equal("Wiwaxia", rownames(read)[4])
  expect_equal("(01)", as.character(read[1, 27]))
})

test_that("TNT trees parsed correctly", {
  trees <- ReadTntTree('test-tnt-tree.tre')
  expect_equal(2, length(trees))
  expect_equal(32, ConsensusWithout(trees, 'Paterimitra')$Nnode)
})
