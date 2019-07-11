library('TreeDist')
trees <- phangorn::allTrees(7, rooted=FALSE, tip.label=letters[1:7])
trees <- lapply(trees, ape::root, 'a', resolve.root=TRUE)
trees <- trees[1:4]
TreeShape7 <- function (tree) {
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  tab <- table(rowSums(TreeSearch::AllDescendantEdges(parent, child, 12)))
  as.integer(sum(tab * (10 ^ ((as.integer(names(tab)) - 1L) / 2))))
}
treeShapes <- vapply(trees, TreeShape7, 1L)
trees <- trees[order(treeShapes)]

# First check that VAI is working!
vai <- VariationOfArborealInfo(trees, trees, normalize=TRUE)

elementStatus <- Quartet::ManyToManyQuartetAgreement(trees)
qd <- 1 - (rowSums(elementStatus[, , c('d', "d", "r1", "r2"), drop = FALSE], dims=2) / 
       rowSums(elementStatus, dims=2))

treeDists <- vapply(trees, function (tr1) vapply(trees, function (tr2) {
  c(phangorn::treedist(tr1, tr2)[c('symmetric.difference', 'path.difference')],
    phangorn::SPR.dist(tr1, tr2))
}, double(3)), matrix(0, nrow=3, ncol=length(trees)))

sevenTipDistances <- list(
 vai = vai,
 vpi = VariationOfPartitionInfo(trees, trees, normalize=TRUE),
 vci = VariationOfClusteringInfo(trees, trees, normalize=TRUE),
 qd = qd,
 nts = 1 - NyeTreeSimilarity(trees, trees, normalize=TRUE),
 msd = MatchingSplitDistance(trees, trees),
 rf = treeDists['symmetric.difference', , ],
 path = treeDists['path.difference', , ],
 spr = treeDists['spr', , ],
 shapes = treeShapes[order(treeShapes)]
)

usethis::use_data(sevenTipDistances, overwrite=TRUE)
