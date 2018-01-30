library(ape)
library(testthat)

## Test cases designed by Thomas Guillerme
context("Morphy: Correct step counting")
test_that("Morphy generates correct lengths", {
  ## Tree
  tree <- ape::read.tree(text = "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));")
  characters <- c("23--1??--032", # 0,  expect score = 5 
                  "1---1111---1", # 1,  expect score = 2
                  "1100----1100", # 2,  expect score = 3
                  "11-------100", # 3,  expect score = 2
                  "----1111---1", # 4,  expect score = 1
                  "01----010101", # 5,  expect score = 5
                  "01---1010101", # 6,  expect score = 5
                  "1??--??--100", # 7,  expect score = 2
                  "21--3??--032", # 8,  expect score = 5
                  "11--1??--111", # 9,  expect score = 2
                  "11--1000001-", # 10, expect score = 2
                  "01------0101", # 11, expect score = 4
                  "110--?---100", # 12, expect score = 3
                  "11--1??--111", # 13, expect score = 2
                  "210--100--21", # 14, expect score = 5
                  "????----1???", # 15, expect score = 0
                  "23--1----032", # 16, expect score = 5
                  "1----1----1-", # 17, expect score = 2
                  "-1-1-1--1-1-", # 18, expect score = 4
                  "23--1??--032", # 19, expect score = 5
                  "--------0101", # 20, expect score = 2
                  "10101-----01", # 21, expect score = 4
                  "011--?--0011", # 22, expect score = 3
                  "110--??--100", # 23, expect score = 3
                  "11--1000001-", # 24, expect score = 2
                  "21--1----012", # 25, expect score = 5
                  "11----111111", # 26, expect score = 1
                  "10101-----01", # 27, expect score = 4
                  "210210------", # 28, expect score = 4
                  "----1111----", # 29, expect score = 0
                  "230--??1--32", # 30, expect score = 5
                  "023--??1--32", # 31, expect score = 5
                  "023-???1--32", # 32, expect score = 4
                  "23--1?1--023", # 33, expect score = 5
                  "----1010----", # 34, expect score = 2
                  "------11---1", # 35, expect score = 1
                  "10----11---1", # 36, expect score = 3
                  "320--??3--21", # 37, expect score = 5
                  "000011110000"  # 38, expect score = 2
                  ) 
  ## Results
  expected_results <- c(5, 2, 3, 2, 1, 5, 5, 2, 5, 2, 2, 4, 3, 2, 5, 0, 5, 2, 4, 5, 2, 4, 3, 3, 2, 5, 1, 4, 4, 0, 5, 5, 4, 5, 2, 1, 3, 5, 2)

  ##plot(tree); nodelabels(12:22); tiplabels(0:11)
  ## Run the tests
  for(test in seq_along(characters)) {
    morphyObj <- SingleCharMorphy(characters[test])
    tree_length <- MorphyTreeLength(tree, morphyObj)
    #if (tree_length != expected_results[test]) cat("Test case", test - 1, characters[test], "unequal: Morphy calcluates",
    #  tree_length, "instead of", expected_results[test],"\n")
    expect_equal(tree_length, expected_results[test])
    morphyObj <- UnloadMorphy(morphyObj)
  }
  ## Test combined matrix
  bigPhy <- StringToPhyDat(paste0(characters, collapse='\n'), tree$tip.label, byTaxon=FALSE)
  expect_identical(characters, PhyToString(bigPhy, byTaxon=FALSE, concatenate=FALSE))
  expect_identical(paste0(collapse='', vapply(characters, substr, start=0, stop=1, character(1))),
                   substr(PhyToString(bigPhy, ps=';', useIndex=TRUE, byTaxon=TRUE, concatenate=TRUE),
                    start=0, stop=length(characters)))
  morphyObj <- PhyDat2Morphy(bigPhy)
  tree_length <- MorphyTreeLength(tree, morphyObj)
  expect_equal(tree_length, sum(expected_results))

  ## Run the bigger tree tests
  bigTree <- read.tree(text = "((1,2),((3,(4,5)),(6,(7,(8,(9,(10,((11,(12,(13,(14,15)))),(16,(17,(18,(19,20))))))))))));")
  bigChars <- c("11111---111---11---1")
  ## Results
  expected_results <- c(3)

  ## Run the tests
  for(test in 1:length(bigChars)) {
    phy <- StringToPhyDat(bigChars[test], bigTree$tip.label)
    # Presently a good test to confirm that PhyDat2Morphy works with single-character phys
    morphyObj <- PhyDat2Morphy(phy)
    tree_length <- MorphyTreeLength(bigTree, morphyObj)
    #if (tree_length != expected_results[test]) cat("Test case", test - 1, bigChars[test], "unequal: Morphy calcluates",
    #  tree_length, "instead of", expected_results[test],"\n")
    expect_equal(tree_length, expected_results[test])
    UnloadMorphy(morphyObj)
  }
})
