context("RWTY")

if (requireNamespace('rwty', quietly = TRUE)) {
  test_that('Fungi MSID/PID differ', {
    data('fungus', package = 'rwty')
    trees <- fungus[[1]]$trees
    PhylogeneticInfoDistance(trees[1:4])
    MatchingSplitInfoDistance(trees[1:4])
  })
}
