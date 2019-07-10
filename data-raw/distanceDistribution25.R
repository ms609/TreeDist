library(TreeSearch)
library(TreeDist)

logDoubleFactorials <- TreeSearch::logDoubleFactorials
cbPalette8 <- Ternary::cbPalette8
suppressWarnings(RNGversion("3.5.0")) # Stopgap until we can require R 3.6.0
set.seed(0)

RandomDistances <- function (nLeaves, repls) {
  RandomTree <- function(nTip) ape::rtree(nTip, br=NULL)
  vapply(seq_len(repls), 
         function (XX) {
           tr1 <- RandomTree(nLeaves) 
           tr2 <- RandomTree(nLeaves) 
           cat('.')
            
           c(VariationOfArborealInfo(tr1, tr2, normalize=TRUE),
             VariationOfPartitionInfo(tr1, tr2, normalize=TRUE),
             VariationOfClusteringInfo(tr1, tr2, normalize=TRUE),
             Quartet::QuartetDivergence(Quartet::QuartetStatus(tr1, tr2), similarity = FALSE),
             1 - NyeTreeSimilarity(tr1, tr2, normalize=TRUE),
             MatchingSplitDistance(tr1, tr2),
             phangorn::treedist(tr1, tr2),
             phangorn::SPR.dist(tr1, tr2)
           )
         },
         c(vai = 0, vpi = 0, vci = 0, qd = 0, nts = 0, 
           msd = 0, rf = 0, path = 0, spr = 0)
  )
}

distanceDistribution25 <- RandomDistances(25, 10000)
usethis::use_data(distanceDistribution25, overwrite=TRUE)
