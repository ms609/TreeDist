#' Information-based generalized Robinson–Foulds distances
#'
#' Calculate tree similarity and distance measures based on the amount of 
#' phylogenetic or clustering information that two trees hold in common, as
#' proposed in Smith (2020).
#' 
#' 
#' [Generalized Robinson–Foulds distances](https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html#generalized-robinson-foulds-distances)
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
#' when they are optimally matched, following Smith (2020).
#' The complementary tree distance measures state how much information is 
#' different in the splits of two trees, under an optimal matching.
#' 
#' # Concepts of information
#'
#' The phylogenetic (Shannon) information content and entropy of a split are
#' defined in 
#' [a separate vignette](https://ms609.github.io/TreeDist/articles/information.html).
#' 
#' Using the mutual (clustering) information (Meil\ifelse{html}{&#259;}{a}
#' 2007, Vinh _et al._ 2010) of two splits to quantify their similarity gives 
#' rise to the Mutual Clustering Information measure (`MutualClusteringInfo()`, 
#' `MutualClusteringInfoSplits()`); the entropy distance 
#' gives the Clustering Information Distance (`ClusteringInfoDistance()`).
#' This approach is optimal in many regards, and is implemented with 
#' normalization in the convenience function `TreeDistance()`.
#' 
#' Using the amount of phylogenetic information common to two splits to measure
#' their similarity gives rise to the Shared Phylogenetic Information similarity
#' measure (`SharedPhylogeneticInfo()`, `SharedPhylogeneticInfoSplits()`).
#' The amount of information distinct to
#' each of a pair of splits provides the complementary Different Phylogenetic
#' Information distance metric (`DifferentPhylogeneticInfo()`).
#' 
#' The Matching Split Information measure (`MatchingSplitInfo()`,
#' `MatchingSplitInfoSplits()`) defines the similarity between a pair of 
#' splits as the phylogenetic information content of the most informative 
#' split that is consistent with both input splits; `MatchingSplitInfoDistance()`
#' is the corresponding measure of tree difference.
#' [(More information here.)](https://ms609.github.io/TreeDist/articles/Generalized-RF.html)
#' 
#' ## Conversion to distances
#' 
#' To convert similarity measures to distances, it is necessary to 
#' subtract the similarity score from a maximum value.  In order to generate
#' distance _metrics_, these functions subtract the similarity twice from the 
#' total information content (SPI, MSI) or entropy (MCI) of all the splits in 
#' both trees (Smith 2020).
#' 
#' ## Normalization
#' 
#' If `normalize = TRUE`, then results will be rescaled such that distance
#' ranges from zero to (in principle) one.
#' The maximum **distance** is the sum of the information content or entropy of
#' each split in each tree; the maximum **similarity** is half this value.
#' (See Vinh _et al._ (2010, table 3) and Smith (2020) for
#' alternative normalization possibilities.)
#' 
#' Note that a distance value of one (= similarity of zero) will seldom be
#' achieved, as even the most different trees exhibit some similarity.
#' It may thus be helpful to rescale the normalized value such that the
#' _expected_ distance between a random pair of trees equals one.  This can
#' be calculated with `ExpectedVariation()`; or see package
#' '[TreeDistData](https://ms609.github.io/TreeDistData/reference/randomTreeDistances.html)'
#' for a compilation of expected values under different metrics for trees with
#' up to 200 leaves.
#' 
#' Alternatively, to scale against the information content or entropy of all 
#' splits in the most or least informative tree, use `normalize = `[`pmax`] or 
#' [`pmin`] respectively.
#' To calculate the relative similarity against a reference tree that is known
#' to be 'correct', use `normalize = ``SplitwiseInfo(trueTree)` (SPI, MSI) or
#' `ClusteringEntropy(trueTree)` (MCI).
#' 
#' @template tree12ListParams
#' 
#' @param normalize If a numeric value is provided, this will be used as a 
#' maximum value against which to rescale results.
#' If `TRUE`, results will be rescaled against a maximum value calculated from
#' the specified tree sizes and topology, as specified in the 'Normalization' 
#' section below.
#' If `FALSE`, results will not be rescaled.
#' 
#' @param diag Logical specifying whether to return similarities along the
#' diagonal, i.e. of each tree with itself.  Applies only if `tree2` is
#' a list identical to `tree1`, or `NULL`.
#' 
#' @param reportMatching Logical specifying whether to return the clade
#' matchings as an attribute of the score.
#'
#' @return If `reportMatching = FALSE`, the functions return a numeric 
#' vector specifying the requested similarities or differences.
#' 
#' If `reportMatching = TRUE`, the functions additionally return details
#' of which clades are matched in the optimal matching, which can be viewed
#' using [`VisualizeMatching()`].
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
#' # Every split in tree1 conflicts with every split in tree3
#' # Pairs of conflicting splits contain clustering, but not phylogenetic, 
#' # information
#' SharedPhylogeneticInfo(tree1, tree3) # = 0
#' MutualClusteringInfo(tree1, tree3) # > 0
#' 
#' # Converting trees to Splits objects can speed up multiple comparisons
#' splits1 <- TreeTools::as.Splits(tree1)
#' splits2 <- TreeTools::as.Splits(tree2)
#' 
#' SharedPhylogeneticInfoSplits(splits1, splits2)
#' MatchingSplitInfoSplits(splits1, splits2)
#' MutualClusteringInfoSplits(splits1, splits2)
#' @template MRS 
#' 
#' @references 
#'  * \insertRef{Day1985}{TreeDist}
#'  
#'  * \insertRef{Mackay2003}{TreeDist}
#'  
#'  * \insertRef{Meila2007}{TreeDist}
#'  
#'  * \insertRef{SmithDist}{TreeDist}
#'  
#'  * \insertRef{Vinh2010}{TreeDist}
#' 
#' @encoding UTF-8
#' @family tree distances
#' @export
TreeDistance <- function (tree1, tree2 = tree1) {
  ClusteringInfoDistance(tree1, tree2, normalize = TRUE, reportMatching = FALSE)
}

#' @rdname TreeDistance
#' @export
SharedPhylogeneticInfo <- function (tree1, tree2 = NULL, normalize = FALSE,
                                    reportMatching = FALSE, diag = TRUE) {
  unnormalized <- CalculateTreeDistance(SharedPhylogeneticInfoSplits, tree1,
                                        tree2, reportMatching = reportMatching)
  
  if (diag && inherits(unnormalized, 'dist')) {
    unnormalized <- as.matrix(unnormalized)
    diag(unnormalized) <- SplitwiseInfo(tree1)
  }
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = SplitwiseInfo, Combine = .PairMean)
}

#' @rdname TreeDistance
#' @export
DifferentPhylogeneticInfo <- function (tree1, tree2 = NULL, normalize = FALSE,
                                       reportMatching = FALSE) {
  spi <- SharedPhylogeneticInfo(tree1, tree2, normalize = FALSE, diag = FALSE,
                                reportMatching = reportMatching)
  treesIndependentInfo <- .MaxValue(tree1, tree2, SplitwiseInfo)
  
  ret <- treesIndependentInfo - spi - spi
  ret <- NormalizeInfo(ret, tree1, tree2, how = normalize, 
                       infoInBoth = treesIndependentInfo,
                       InfoInTree = SplitwiseInfo, Combine = '+')
  
  ret[ret < .Machine$double.eps ^ 0.5] <- 0 # Catch floating point inaccuracy
  attributes(ret) <- attributes(spi)
  
  # Return:
  ret
}

#' @rdname TreeDistance
#' @export
PhylogeneticInfoDistance <- DifferentPhylogeneticInfo

#' @rdname TreeDistance
#' @aliases ClusteringInfoDist
#' @export
ClusteringInfoDistance <- function (tree1, tree2 = NULL, normalize = FALSE,
                                       reportMatching = FALSE) {
  mci <- MutualClusteringInfo(tree1, tree2, normalize = FALSE, diag = FALSE,
                              reportMatching = reportMatching)
  treesIndependentInfo <- .MaxValue(tree1, tree2, ClusteringEntropy)
  
  ret <- treesIndependentInfo - mci - mci
  ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                       infoInBoth = treesIndependentInfo,
                       InfoInTree = ClusteringEntropy, Combine = '+')
  
  ret[ret < .Machine$double.eps^0.5] <- 0 # In case of floating point inaccuracy
  attributes(ret) <- attributes(mci)
  
  # Return:
  ret
}

# TODO check that Internalizing this hasn't hidden TreeDistance from the index.
#' @export
ClusteringInfoDist <- ClusteringInfoDistance

#' @rdname TreeDistance
#' @param samples Integer specifying how many samplings to obtain; 
#' accuracy of estimate increases with `sqrt(samples)`.
#' @importFrom stats sd
#' @importFrom TreeTools as.Splits
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
      MatchingSplitInfoSplits(splits1, resampled2),
      MutualClusteringInfoSplits(splits1, resampled2)
      )
  }, c(SharedPhylogeneticInfo = 0, MatchingSplitInfo = 0,
       MutualClusteringInfo = 0))
  
  mut <- cbind(Estimate = rowMeans(mutualEstimates),
               sd = apply(mutualEstimates, 1, sd), n = samples)
  
  ret <- rbind(
    mut,
    DifferentPhylogeneticInfo = c(info1 + info2 - mut[1, 1] - mut[1, 1],
                                    mut[1, 2] * 2, samples),
    MatchingSplitInfoDistance = c(info1 + info2 - mut[2, 1] - mut[2, 1],
                                     mut[2, 2] * 2, samples),
    ClusteringInfoDistance = c(ClusteringEntropy(tree1) + 
                                 ClusteringEntropy(tree2) - 
                                 mut[3, 1] - mut[3, 1],
                               mut[3, 2] * 2, samples)
  )
  
  # Return:
  cbind(Estimate = ret[, 1],
        'Std. Err.' = ret[, 'sd'] / sqrt(ret[, 'n']), 
        ret[, 2:3])
}

#' @rdname TreeDistance
#' @aliases MutualClusteringInformation
#' @export
MutualClusteringInfo <- function (tree1, tree2 = NULL, normalize = FALSE,
                                  reportMatching = FALSE, diag = TRUE) {
  unnormalized <- CalculateTreeDistance(MutualClusteringInfoSplits, tree1, tree2,
                                        reportMatching)
  if (diag && is.null(tree2)) {
    unnormalized <- as.matrix(unnormalized)
    diag(unnormalized) <- ClusteringEntropy(tree1)
  }
  NormalizeInfo(unnormalized, tree1, tree2, ClusteringEntropy, 
                how = normalize, Combine = .PairMean)
}

#' @export
MutualClusteringInformation <- MutualClusteringInfo

#' @rdname TreeDistance
#' @template splits12params
#' @param nTip (Optional) Integer specifying the number of leaves in each split.
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

#' Mean of two numbers
#' 
#' Used for normalization and range calculation
#' 
#' @keywords internal
.PairMean <- function (x, y) (x + y) / 2L
