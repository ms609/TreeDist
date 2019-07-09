library(TreeSearch)
cbPalette8 <- Ternary::cbPalette8
nLeaves <- c(seq(5, 50, by=5), seq(60, 100, by=10), seq(150, 200, by=25))

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
                            1 - NyeTreeSimilarity(tr1, tr2, normalize=TRUE)
                            )
                        },
                        double(5))
    cbind(rowMeans(distances), apply(distances, 1, sd))
    }, matrix(0, ncol=2, nrow=5, dimnames=list(c('vai', 'vpi', 'vci', 'qd', 'nts'), c('mean', 'sd')))
  )
  dimnames(ret)[[3]] <- nLeaves
  ret
}

AddLine <- function (dat, col) {
  nLeaves <- colnames(dat)
  points(nLeaves, dat['mean', ], pch=3, col=col)
  lines(nLeaves, dat['mean', ] + dat['sd', ], lty=3, col=col)
  lines(nLeaves, dat['mean', ] - dat['sd', ], lty=3, col=col)
}


dists <- RandomDistances(nLeaves, repls=100L)

plot(nLeaves, dists['vai', 'mean', ], pch='.', ylim=c(0.3, 1), col='white')
AddLine(dists['vai', , ], col=cbPalette8[1])
AddLine(dists['vpi', , ], col=cbPalette8[2])
AddLine(dists['vci', , ], col=cbPalette8[3])
AddLine(dists['qd', , ], col=cbPalette8[4])
AddLine(dists['nts', , ], col=cbPalette8[5])

legend('bottomright', bty='n', pch=3, legend=dimnames(dists)[[1]], col=cbPalette8[1:5])

randomTreeDistances <- dists
usethis::use_data(randomTreeDistances, overwrite=TRUE)
