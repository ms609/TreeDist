library(TreeDist)
suppressWarnings(RNGversion("3.5.0")) # Stopgap until we can require R 3.6.0
set.seed(0)

RandomDistances <- function (nLeaves, repls) {
  RandomTree <- function(nTip) ape::rtree(nTip, br=NULL)
  ret <- vapply(nLeaves, function (n) {
    cat('\n', n, 'Leaves ')
    distances <- vapply(seq_len(repls), 
                        function (XX) {
                          tr1 <- RandomTree(n) 
                          tr2 <- RandomTree(n) 
                          cat('.')
                          
                          c(VariationOfArborealInfo(tr1, tr2, normalize=TRUE),
                            VariationOfPartitionInfo(tr1, tr2, normalize=TRUE),
                            VariationOfClusteringInfo(tr1, tr2, normalize=TRUE),
                            Quartet::QuartetDivergence(Quartet::QuartetStatus(tr1, tr2), similarity = FALSE),
                            1 - NyeTreeSimilarity(tr1, tr2, normalize=TRUE),
                            MatchingSplitDistance(tr1, tr2, normalize=FALSE),
                            phangorn::treedist(tr1, tr2) / c(n + n - 6L, 1), # No norm for path
                            phangorn::SPR.dist(tr1, tr2) / (n / 2) # crude normalization!
                            )
                        },
                        double(9))
    cbind(rowMeans(distances), apply(distances, 1, sd))
    }, matrix(0, ncol=2, nrow=9, dimnames=list(c('vai', 'vpi', 'vci', 'qd', 'nts', 
                                                 'msd', 'rf', 'path', 'spr'), c('mean', 'sd')))
  )
  dimnames(ret)[[3]] <- nLeaves
  ret
}

nLeaves <- c(seq(5, 50, by=5), seq(60, 90, by=10), seq(100, 200, by=20))
randomTreeDistances <- RandomDistances(nLeaves, repls=3L)
#usethis::use_data(randomTreeDistances, compress='gzip', overwrite=TRUE)
