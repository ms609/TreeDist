context('Tree differences')

test_that("Split combatibility is correctly established", {
  expect_true(SplitsCompatible(as.logical (c(0,0,1,1,0)), as.logical(c(0,0,1,1,0))))
  expect_true(SplitsCompatible(as.logical (c(0,0,1,1,0)), !as.logical(c(0,0,1,1,0))))
  expect_true(SplitsCompatible(as.logical (c(0,0,1,1,0)), as.logical(c(1,0,1,1,0))))
  expect_true(SplitsCompatible(!as.logical(c(0,0,1,1,0)), as.logical(c(1,0,1,1,0))))
  expect_false(SplitsCompatible(as.logical(c(0,0,1,1,0)), as.logical(c(1,1,0,1,0))))
})

test_that('Tree differences are correctly calculated', {
  # Labels in different order to confound Tree2Splits
  treeSym8 <- ape::read.tree(text='((e, (f, (g, h))), (((a, b), c), d));')
  treeBal8 <- ape::read.tree(text='(((e, f), (g, h)), ((a, b), (c, d)));')
  treeOpp8 <- ape::read.tree(text='(((a, f), (c, h)), ((g, b), (e, d)));')
  treeBadLabel8 <- ape::read.tree(text='((a, b, c, D), (e, f, g, h));')
  
  treeCat8 <- ape::read.tree(text='((((h, g), f), e), (d, (c, (b, a))));')
  treeTac8 <- ape::read.tree(text='((((e, c), g), a), (h, (b, (d, f))));')
  
  treeAb.Cdefgh <- ape::read.tree(text='((a, b), (c, d, e, f, g, h));')
  treeAbc.Defgh <- ape::read.tree(text='((a, b, c), (d, e, f, g, h));')
  treeAcd.Befgh <- ape::read.tree(text='((a, c, d), (b, e, f, g, h));')
  treeAbcd.Efgh <- ape::read.tree(text='((a, b, c, d), (e, f, g, h));')
  treeTwoSplits <- ape::read.tree(text="(((a, b), c, d), (e, f, g, h));")

  # Labels differ
  expect_error(MutualArborealInfo(treeSym8, treeBadLabel8))
  expect_error(MutualPartitionInfo(treeSym8, treeBadLabel8))
  expect_error(VariationOfArborealInfo(treeSym8, treeBadLabel8))
  expect_error(VariationOfPartitionInfo(treeSym8, treeBadLabel8))
  expect_error(MutualClusterInfo(treeSym8, treeBadLabel8))
  expect_error(NyeTreeDistance(treeSym8, treeBadLabel8))
  expect_error(MatchingSplitDistance(treeSym8, treeBadLabel8))

  expect_equal(22.53747, MutualArborealInfo(treeSym8, treeSym8), tolerance=1e-05)
  expect_equal(13.75284, MutualArborealInfo(treeSym8, treeBal8), tolerance=1e-05)
  expect_equal(VariationOfArborealInfo(treeSym8, treeAcd.Befgh),
               VariationOfArborealInfo(treeAcd.Befgh, treeSym8), tolerance=1e-05)
  expect_equal(22.53747 + MutualArborealInfo(treeAcd.Befgh, treeAcd.Befgh) - 
                 (2 * MutualArborealInfo(treeSym8, treeAcd.Befgh)), 
               VariationOfArborealInfo(treeSym8, treeAcd.Befgh), tolerance=1e-05)
  expect_equal(-log2(945/10395), MutualArborealInfo(treeSym8, treeAb.Cdefgh))
  expect_equal(22.53747 + MutualArborealInfo(treeBal8, treeBal8) - 13.75284 - 13.75284, 
               VariationOfArborealInfo(treeSym8, treeBal8), tolerance=1e-05)
  expect_equal(-log2(945/10395), MutualArborealInfo(treeSym8, treeAb.Cdefgh))
  expect_equal(-log2(315/10395), MutualArborealInfo(treeSym8, treeAbc.Defgh))
  expect_equal(0, VariationOfArborealInfo(treeSym8, treeSym8))
  expect_equal(PartitionInfo(treeSym8) - PartitionInfo(treeAcd.Befgh),
               VariationOfArborealInfo(treeSym8, treeAbc.Defgh))
  
  BinaryToSplit <- function (binary) matrix(as.logical(binary))
  expect_equal(MutualPartitionInfoSplits(
    BinaryToSplit(c(1, 1, 0, 0, 0, 0, 0, 0)),
    BinaryToSplit(c(0, 0, 1, 1, 0, 0, 0, 0))
    ), MutualPartitionInfoSplits(
    BinaryToSplit(c(0, 0, 0, 0, 0, 0, 1, 1)),
    BinaryToSplit(c(0, 0, 1, 1, 0, 0, 0, 0))
    ))
  
  MutualPartitionInfoSplits(BinaryToSplit(c(1, 1, 1, 1, 0, 0, 0, 0)),
                         BinaryToSplit(c(1, 0, 1, 0, 1, 0, 1, 0)))
  expect_equal(MutualArborealInfo(treeSym8, treeSym8),
               MutualPartitionInfo(treeSym8, treeSym8), tolerance=1e-05)
  expect_equal(MutualPartitionInfo(treeAb.Cdefgh, treeAbc.Defgh),
               MutualPartitionInfo(treeAbc.Defgh, treeAb.Cdefgh))
  expect_equal(MutualPartitionInfo(treeAbcd.Efgh, treeAb.Cdefgh),
               MutualPartitionInfo(treeAb.Cdefgh, treeAbcd.Efgh))
  expect_equal(-(LogTreesMatchingSplit(2, 5) - LnUnrooted.int(7)) / log(2), 
               MutualPartitionInfo(treeAb.Cdefgh, treeAbc.Defgh))
  expect_true(MutualPartitionInfo(treeSym8, treeBal8) > MutualPartitionInfo(treeSym8, treeOpp8))
  expect_equal(0, VariationOfPartitionInfo(treeSym8, treeSym8))
  
  
  expect_equal(5L, NyeTreeSimilarity(treeSym8, treeSym8))
  expect_equal(2, 3 * NyeTreeSimilarity(treeAb.Cdefgh, treeAbc.Defgh))
  expect_equal(3.8, NyeTreeSimilarity(treeSym8, treeBal8))
  expect_true(NyeTreeSimilarity(treeSym8, treeBal8) > NyeTreeSimilarity(treeSym8, treeOpp8))
  
  expect_equal(0L, MatchingSplitDistance(treeSym8, treeSym8))
  expect_equal(1L, MatchingSplitDistance(treeAb.Cdefgh, treeAbc.Defgh))
  expect_equal(2L, MatchingSplitDistance(treeAb.Cdefgh, treeAbcd.Efgh))
  
  # Invariant to tree description order
  sq_pectinate <- ape::read.tree(text='((((((1, 2), 3), 4), 5), 6), (7, (8, (9, (10, 11)))));')
  shuffle1 <- ape::read.tree(text='(((((1, 5), 2), 6), (3, 4)), ((8, (7, 9)), (10, 11)));')
  shuffle2 <- ape::read.tree(text='(((8, (7, 9)), (10, 11)), ((((1, 5), 2), 6), (3, 4)));')
  expect_equal(MatchingSplitDistance(shuffle1, sq_pectinate),
               MatchingSplitDistance(sq_pectinate, shuffle1))
  expect_equal(0L, MatchingSplitDistance(shuffle1, shuffle2))
  expect_equal(MatchingSplitDistance(shuffle1, sq_pectinate),
               MatchingSplitDistance(shuffle2, sq_pectinate))
  
  expect_true(NyeTreeSimilarity(treeSym8, treeBal8) > NyeTreeSimilarity(treeSym8, treeOpp8))
  
  # Test symmetry of small vs large splits
  expect_equal(MutualArborealInfo(treeSym8, treeAbc.Defgh),
               MutualArborealInfo(treeAbc.Defgh, treeSym8))
  expect_equal(-log2(225/10395), MutualArborealInfo(treeSym8, treeAbcd.Efgh))
  expect_equal(-log2(225/10395) - log2(945/10395),
               MutualArborealInfo(treeSym8, treeTwoSplits))
  expect_equal(SplitMutualInformation(8, 4, 3),
               MutualArborealInfo(treeTwoSplits, treeAbc.Defgh))
  
  expect_equal(MutualArborealInfo(treeSym8, list(treeSym8, treeBal8)), 
               MutualArborealInfo(list(treeSym8, treeBal8), treeSym8))
  expect_equal(matrix(c(MutualArborealInfo(treeSym8, treeSym8),
                        MutualArborealInfo(treeBal8, treeSym8),
                        MutualArborealInfo(treeSym8, treeAbc.Defgh),
                        MutualArborealInfo(treeBal8, treeAbc.Defgh),
                        MutualArborealInfo(treeSym8, treeAbcd.Efgh), 
                        MutualArborealInfo(treeBal8, treeAbcd.Efgh)),
                        3L, 2L, byrow=TRUE,
                        dimnames=list(c('sym', 'abc', 'abcd'), c('sym', 'bal'))), 
    MutualArborealInfo(list(sym=treeSym8, bal=treeBal8), 
               list(sym=treeSym8, abc=treeAbc.Defgh, abcd=treeAbcd.Efgh)))
  
})
