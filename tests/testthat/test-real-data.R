#context("RWTY")
#
#if (requireNamespace('rwty', quietly = TRUE)) {
#  test_that('Fungi MSID/PID differ', {
#    data('fungus', package = 'rwty')
#    trees <- fungus[[1]]$trees
#    
#    expect_equal(PhylogeneticInfoDistance(trees[[1]], trees[[2]]),
#                 PhylogeneticInfoDistance(trees[1:2])[1])
#    expect_equal(MatchingSplitInfoDistance(trees[[1]], trees[[2]]),
#                 MatchingSplitInfoDistance(trees[1:2])[1])
#    expect_false(identical(PhylogeneticInfoDistance(trees[1:4]),
#                           MatchingSplitInfoDistance(trees[1:4])))
#  })
#  
#  if (requireNamespace('whalehead', quietly = TRUE)) {
#    test_that('Whalehead MSID/PID differ', {
#      whaleRoot <- system.file('Data', package = 'whalehead', mustWork = TRUE)
#      whaleRoot <- '../ProjectWhalehead/Data'
#      treeFiles <- list.files(path = paste0(whaleRoot, "/CombinedTrees"),
#                              pattern = "*.nex",
#                              full.names = TRUE, recursive = FALSE)
#      trees <- unname(rwty::load.trees(treeFiles[2], trim = 100, skip = 0)$trees)
#      
#      expect_equal(PhylogeneticInfoDistance(trees[[1]], trees[[2]]),
#                   PhylogeneticInfoDistance(trees[1:2])[1])
#      expect_equal(MatchingSplitInfoDistance(trees[[1]], trees[[2]]),
#                   MatchingSplitInfoDistance(trees[1:2])[1])
#      expect_false(identical(PhylogeneticInfoDistance(trees[1:4]),
#                             MatchingSplitInfoDistance(trees[1:4])))
#    })
#  }
#}
#