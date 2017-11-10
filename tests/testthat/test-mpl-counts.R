library(ape)
library(testthat)

## Test cases designed by Thomas Guillerme
context("Correct step counting")
test_that("Morphy generates correct lengths", {
  ## Tree
  tree <- read.tree(text = "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));")
  characters <- c("23--1??--032", # 0,  expect score = 5
                  "1---1111---1", # 1,  expect score = 2
                  "1100----1100", # 2,  expect score = 3
                  "11-------100", # 3,  expect score = 2
                  "----1111---1", # 4,  expect score = 1
                  "01----010101", # 5,  expect score = 5
                  "01---1010101", # 6,  expect score = 5
                  "1??--??--100", # 7,  expect score = 2
                  "21--3??--032", # 8,  expect score = 2
                  "11--1??--111", # 9,  expect score = 2
                  "11--1000001-", # 10, expect score = 4
                  "01------0101", # 11
                  "110--?---100", # 12
                  "11--1??--111", # 13
                  "210--100--21", # 14
                  "????----1???", # 15
                  "23--1----032", # 16
                  "1----1----1-", # 17
                  "-1-1-1--1-1-", # 18
                  "23--1??--032", # 19
                  "--------0101", # 20
                  "10101-----01", # 21
                  "011--?--0011", # 22
                  "110--??--100", # 23
                  "11--1000001-", # 24
                  "21--1----012", # 25
                  "11----111111", # 26
                  "10101-----01", # 27
                  "210210------", # 28
                  "----1111----", # 29
                  "230--??1--32", # 30
                  "023--??1--32", # 31
                  "023-???1--32", # 32
                  "23--1?1--023", # 33
                  "----1010----", # 34
                  "------11---1", # 35
                  "10----11---1", # 36
                  "320--??3--21", # 37
                  "000011110000"  # 38
                  ) 
  ## Results
  expected_results <- c(5, 2, 3, 2, 1, 5, 5, 2, 5, 2, 2, 4, 3, 2, 5, 0, 5, 2, 4, 5, 2, 4, 3, 3, 2, 5, 1, 4, 4, 0, 5, 5, 4, 5, 2, 1, 3, 5, 2)

  ##plot(tree); nodelabels(12:22); tiplabels(0:11)
  ## Run the tests
  for(test in seq_along(characters)) {
    phy <- StringToPhyDat(characters[test], tree$tip.label)
    morphyObj <- PhyDat2Morphy(phy)
    tree_length <- MorphyLength(tree, morphyObj)
    #if (tree_length != expected_results[test]) cat("Test case", test - 1, characters[test], "unequal: Morphy calcluates",
    #  tree_length, "instead of", expected_results[test],"\n")
    expect_equal(tree_length, expected_results[test])
    morphyObj <- UnloadMorphy(morphyObj)
  }
  ## Test combined matrix
  bigPhy <- StringToPhyDat(paste0(characters, collapse='\n'), tree$tip.label, byTaxon=FALSE)
  morphyObj <- PhyDat2Morphy(bigPhy)
  tree_length <- MorphyLength(tree, morphyObj)
  expect_equal(tree_length, sum(expected_results))

  ## Run the bigger tree tests
  tree <- read.tree(text = "((1,2),((3,(4,5)),(6,(7,(8,(9,(10,((11,(12,(13,(14,15)))),(16,(17,(18,(19,20))))))))))));")
  characters <- c("11111---111---11---1" # 1
                  )
  ## Results
  expected_results <- c(3)

  ## Run the tests
  for(test in 1:length(characters)) {
    phy <- StringToPhyDat(characters[test], tree$tip.label)
    morphyObj <- PhyDat2Morphy(phy)
    tree_length <- MorphyLength(tree, morphyObj)
    #if (tree_length != expected_results[test]) cat("Test case", test - 1, characters[test], "unequal: Morphy calcluates",
    #  tree_length, "instead of", expected_results[test],"\n")
    expect_equal(tree_length, expected_results[test])
    UnloadMorphy(morphyObj)
  }
})
