#' Information-based generalized Robinson-Foulds distances
#'
#' Tree similarity and distance measures that measure the amount of 
#' phylogenetic or clustering information that two trees hold in common.
#' 
#' 
#' [Generalized Robinson-Foulds distances](https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html#generalized-robinson-foulds-distances).
#' calculate tree similarity by finding an
#' optimal matching that the similarity between a split on one tree
#' and its pair on a second, considering all possible ways to pair splits 
#' between trees (including leaving a split unpaired).
#' 
#' The methods implemented here use the concepts of 
#' [entropy and information](https://ms609.github.io/TreeDist/articles/information.html)
#' (MacKay 2003) to assign a similarity score between each pair of splits.
#' 
#' The returned tree similarity measures state the amount of information, 
#' in bits, that the splits in two trees hold in common 
#' when they are optimally matched, following Smith (forthcoming).
#' The complementary tree distance measures state how much information is 
#' different in the splits of two trees, under an optimal matching.
#' 
#' @section Concepts of information:
#'
#' The phylogenetic (Shannon) information content and entropy of a split are
#' defined in 
#' [a separate vignette](https://ms609.github.io/TreeDist/articles/information.html).
#' 
#' Using the mutual (clustering) information (Meila 2007, Vinh _et al._ 2010) of
#' two splits to quantify their similarity gives rise to the Mutual Clustering 
#' Information measure (`MutualClusteringInfo`); the entropy distance 
#' gives the Clustering Distance `ClusteringInfoDistance`).
#' This approach is optimal in many regards, and is implemented, normalized
#' against the total information present, in the convenience function `TreeDistance`.
#' 
#' Using the amount of phylogenetic information common to two splits to measure
#' their similarity gives rise to the Shared Phylogenetic Information similarity
#' measure (`SharedPhylogeneticInfo`).  The amount of information distinct to
#' each of a pair of splits provides the complementary Different Phylogenetic
#'  Information distance metric (`DifferentPhylogeneticInfo`).
#' 
#' The `MatchingSplitInfo` measure defines the similarity between a pair of 
#' splits as the phylogenetic information content of the most informative 
#' split that is consistent with both input splits; `MatchingSplitInfoDistance`
#' is the corresponding measure of tree difference.
#' [(More information here.)](https://ms609.github.io/TreeDist/articles/Generalized-RF.html)
#' 
#' 
#' @section Normalization:
#' 
#' If `normalize = TRUE`, then results will be rescaled from zero to a nominal
#' maximum value, calculated thus:
#' 
#' * `MutualClusteringInfo`, `ClusteringInfoDistance`: The sum of the entropy of
#' each split in each of the two trees.  See Vinh _et al._ (2010, table 3) for
#' alternative normalization variants.
#' 
#' * `SharedPhylogeneticInfo`, `MatchingSplitInfo`:
#'  The sum of the information content of all splits in the least informative
#'  tree. To scale against the information content of all splits in the most
#'  informative tree, use `normalize = pmax`.
#' 
#' * `DifferentPhylogeneticInfo`, `MatchingSplitInfoDistance`: The sum of the
#' phylogenetic information content of each split in each of the two trees.
#' 
#' @template tree12listparams
#' 
#' @param normalize If a numeric value is provided, this will be used as a 
#' maximum value against which to rescale results.
#' If `TRUE`, results will be rescaled against a maximum value calculated from
#' the specified tree sizes and topology, as specified in the 'Normalization' 
#' section below.
#' If `FALSE`, results will not be rescaled.
#' 
#' @param reportMatching Logical specifying whether to return the clade
#' matchings as an attribute of the score.
#'
#' @return If `reportMatching = FALSE`, the functions return a numeric 
#' vector specifying the requested similarities or differences.
#' 
#' If `reportMatching = TRUE`, the functions additionally return details
#' of which clades are matched in the optimal matching, which can be viewed
#' using [`VisualizeMatching`].
#'  
#' @examples 
#' tree1 <- ape::read.tree(text='((((a, b), c), d), (e, (f, (g, h))));')
#' tree2 <- ape::read.tree(text='(((a, b), (c, d)), ((e, f), (g, h)));')
#' tree3 <- ape::read.tree(text='((((h, b), c), d), (e, (f, (g, a))));')
#' 
#' # Best possible score is obtained by matching a tree with itself
#' DifferentPhylogeneticInfo(tree1, tree1) # 0, by definition
#' SharedPhylogeneticInfo(tree1, tree1)
#' SplitwiseInfo(tree1) # Maximum shared phylogenetic information
#' 
#' # Best possible score is a function of tree shape; the splits within
#' # balanced trees are more independent and thus contain less information
#' SplitwiseInfo(tree2)
#' 
#' # How similar are two trees?
#' SharedPhylogeneticInfo(tree1, tree2) # Amount of phylogenetic information in common
#' VisualizeMatching(SharedPhylogeneticInfo, tree1, tree2) # Which clades are matched?
#' 
#' DifferentPhylogeneticInfo(tree1, tree2) # Distance measure
#' DifferentPhylogeneticInfo(tree2, tree1) # The metric is symmetric
#'
#' # Are they more similar than two trees of this shape would be by chance?
#' ExpectedVariation(tree1, tree2, sample=12)['DifferentPhylogeneticInfo', 'Estimate']
#' 
#' # Every split in tree1 is contradicted by every split in tree3
#' # Non-arboreal matches contain clustering, but not phylogenetic, information
#' SharedPhylogeneticInfo(tree1, tree3) # = 0
#' MutualClusteringInfo(tree1, tree3) # > 0
#' 
#' 
#' @references 
#'  * \insertRef{Mackay2003}{TreeDist}
#'  
#'  * \insertRef{Meila2007}{TreeDist}
#'  
#'  * \insertRef{SmithDist}{TreeDist}
#'  
#'  * \insertRef{Vinh2010}{TreeDist}
#' 
#' @template MRS
#' 
#' @family tree distances
#' @export
TreeDistance <- function (tree1, tree2 = tree1) {
  MutualClusteringInfo(tree1, tree2, normalize = TRUE, reportMatching = FALSE)
}

#' @describeIn TreeDistance Shared phylogenetic information between splits of
#'  two trees.
#' @export
SharedPhylogeneticInfo <- function (tree1, tree2 = tree1, normalize = FALSE,
                                    reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(SharedPhylogeneticInfoSplits, tree1,
                                        tree2, reportMatching=reportMatching)
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = SplitwiseInfo, Combine = pmin)
}

#' @describeIn TreeDistance Phylogenetic information difference between splits
#' of two trees.
#' @export
DifferentPhylogeneticInfo <- function (tree1, tree2 = tree1, 
                                         normalize = FALSE,
                                         reportMatching = FALSE) {
  spi <- SharedPhylogeneticInfo(tree1, tree2, normalize = FALSE,
                                reportMatching = reportMatching)
  treesIndependentInfo <- outer(SplitwiseInfo(tree1), SplitwiseInfo(tree2), '+')
  
  ret <- treesIndependentInfo - spi - spi
  ret <- NormalizeInfo(ret, tree1, tree2, how=normalize, 
                       infoInBoth = treesIndependentInfo,
                       InfoInTree = SplitwiseInfo, Combine = '+')
  
  ret[ret < 1e-13] <- 0 # In case of floating point inaccuracy
  attributes(ret) <- attributes(spi)
  
  # Return:
  ret
}

#' @describeIn TreeDistance Variation of clustering information distance between
#' splits of two trees.
#' @aliases ClusteringInfoDist
#' @export
ClusteringInfoDistance <- function (tree1, tree2 = tree1, normalize = FALSE,
                                       reportMatching = FALSE) {
  mci <- MutualClusteringInfo(tree1, tree2, normalize = FALSE, 
                              reportMatching = reportMatching)
  treesIndependentInfo <- outer(ClusteringEntropy(tree1), ClusteringEntropy(tree2), '+')
  
  ret <- treesIndependentInfo - mci - mci
  ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                       infoInBoth = treesIndependentInfo,
                       InfoInTree = ClusteringEntropy, Combine = '+')
  
  ret[ret < 1e-13] <- 0 # In case of floating point inaccuracy
  attributes(ret) <- attributes(mci)
  
  # Return:
  ret
}

# TODO check that Internalizing this hasn't hidden TreeDistance from the index.
#' @export
ClusteringInfoDist <- ClusteringInfoDistance

#' @describeIn TreeDistance Estimate expected Shared and Different
#' Phylogenetic Information for a pair of trees of a given topology.
#' @param samples Integer specifying how many samplings to obtain; 
#' accuracy of estimate increases with `sqrt(samples)`.
#' @importFrom stats sd
#' @importFrom TreeTools as.Splits .DecodeBinary
#' @export
ExpectedVariation <- function (tree1, tree2, samples = 1e+4) {
  info1 <- SplitwiseInfo(tree1)
  info2 <- SplitwiseInfo(tree2)
  splits1 <- as.Splits(tree1)
  tipLabels <- attr(splits1, 'tip.label')
  nTip <- attr(splits1, 'nTip')
  splits2 <- as.Splits(tree2, tipLabels)
  
  mutualEstimates <- vapply(seq_len(samples), function (x) {
    resampled2 <- as.Splits(splits2, sample(tipLabels))
    
    c(SharedPhylogeneticInfoSplits(splits1, resampled2),
      MatchingSplitInfoSplits(splits1, resampled2))
  }, c(SharedPhylogeneticInfo = 0, MatchingSplitInfo = 0))
  
  mut <- cbind(Estimate = rowMeans(mutualEstimates),
               sd = apply(mutualEstimates, 1, sd), n = samples)
  
  ret <- rbind(
    mut,
    DifferentPhylogeneticInfo = c(info1 + info2 - mut[1, 1] - mut[1, 1],
                                    mut[1, 2] * 2, samples),
    MatchingSplitInfoDistance = c(info1 + info2 - mut[2, 1] - mut[2, 1],
                                     mut[2, 2] * 2, samples)
  )
  
  # Return:
  cbind(Estimate = ret[, 1],
        'Std. Err.' = ret[, 'sd'] / sqrt(ret[, 'n']), 
        ret[, 2:3])
}

#' @describeIn TreeDistance Mutual Clustering Information of two trees.
#' @aliases MutualClusteringInformation
#' @export
MutualClusteringInfo <- function (tree1, tree2 = tree1, normalize = FALSE,
                                  reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(MutualClusteringInfoSplits, tree1, tree2,
                                        reportMatching)
  NormalizeInfo(unnormalized, tree1, tree2, ClusteringEntropy, 
                how = normalize, Combine = pmin)
}

#' @export
MutualClusteringInformation <- MutualClusteringInfo

#' @rdname TreeDistance
#' @template splits12params
#' @template nTipParam
#' @export
SharedPhylogeneticInfoSplits <- function (splits1, splits2,
                                          nTip = attr(splits1, 'nTip'),
                                          reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, nTip, cpp_shared_phylo,
                maximize = TRUE, reportMatching = reportMatching)
}

#' @rdname TreeDistance
#' @export
MutualClusteringInfoSplits <- function (splits1, splits2,
                                        nTip = attr(splits1, 'nTip'),
                                        reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, nTip, cpp_mutual_clustering,
                maximize = TRUE, reportMatching = reportMatching)
}
