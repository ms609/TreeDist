#' Information-based generalized Robinson&ndash;Foulds distances
#'
#' Calculate tree similarity and distance measures based on the amount of 
#' phylogenetic or clustering information that two trees hold in common, as
#' proposed in Smith (2020).
#' 
#' 
#' [Generalized Robinson&ndash;Foulds distances](https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html#generalized-robinson-foulds-distances)
#' calculate tree similarity by finding an
#' optimal matching that the similarity between a split on one tree
#' and its pair on a second, considering all possible ways to pair splits 
#' between trees (including leaving a split unpaired).
#' 
#' The methods implemented here use the concepts of 
#' [entropy and information](https://ms609.github.io/TreeDist/articles/information.html)
#' \insertCite{Mackay2003}{TreeDist} to assign a similarity score between each
#' pair of splits.
#' 
#' The returned tree similarity measures state the amount of information, 
#' in bits, that the splits in two trees hold in common 
#' when they are optimally matched, following 
#' \insertCite{SmithDist;textual}{TreeDist}.
#' The complementary tree distance measures state how much information is 
#' different in the splits of two trees, under an optimal matching.
#' Where trees contain different tips, tips present in one tree but not the
#' other are removed before each comparison (as by definition, the trees neither
#' hold information in common nor differ regarding these tips).
#' 
#' # Concepts of information
#'
#' The phylogenetic (Shannon) information content and entropy of a split are
#' defined in 
#' [a separate vignette](https://ms609.github.io/TreeDist/articles/information.html).
#' 
#' Using the mutual (clustering) information
#' \insertCite{Meila2007,Vinh2010}{TreeDist} of two splits to quantify their
#' similarity gives rise to the Mutual Clustering Information measure
#' (`MutualClusteringInfo()`, `MutualClusteringInfoSplits()`);
#' the entropy distance gives the Clustering Information Distance
#' (`ClusteringInfoDistance()`).
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
#' ([More information here](
#' https://ms609.github.io/TreeDist/articles/Generalized-RF.html).)
#' 
#' # Conversion to distances
#' 
#' To convert similarity measures to distances, it is necessary to 
#' subtract the similarity score from a maximum value.  In order to generate
#' distance _metrics_, these functions subtract the similarity twice from the 
#' total information content (SPI, MSI) or entropy (MCI) of all the splits in 
#' both trees \insertCite{SmithDist}{TreeDist}.
#' 
#' # Normalization
#' 
#' If `normalize = TRUE`, then results will be rescaled such that distance
#' ranges from zero to (in principle) one.
#' The maximum **distance** is the sum of the information content or entropy of
#' each split in each tree; the maximum **similarity** is half this value.
#' (See Vinh _et al._ (2010, table 3) and 
#' \insertCite{SmithDist;textual}{TreeDist} for
#' alternative normalization possibilities.)
#' 
#' Note that a distance value of one (= similarity of zero) will seldom be
#' achieved, as even the most different trees exhibit some similarity.
#' It may thus be helpful to rescale the normalized value such that the
#' _expected_ distance between a random pair of trees equals one.  This can
#' be calculated with `ExpectedVariation()`; or see package
#' '[TreeDistData](
#' https://ms609.github.io/TreeDistData/reference/randomTreeDistances.html)'
#' for a compilation of expected values under different metrics for trees with
#' up to 200 leaves.
#' 
#' Alternatively, use `normalize = `[`pmax`] or [`pmin`] to scale against the
#' information content or entropy of all splits in the most (`pmax`) or
#' least (`pmin`) informative tree in each pair.
#' To calculate the relative similarity against a reference tree that is known
#' to be "correct", use `normalize = SplitwiseInfo(trueTree)` (SPI, MSI) or
#' `ClusteringEntropy(trueTree)` (MCI).
#' For worked examples, see the internal function [`NormalizeInfo()`], which is
#' called from distance functions with the parameter `how = normalize`.
#' .
#' 
#'
#' # Distances between large trees
#' 
#' To balance memory demands and runtime with flexibility, these functions are
#' implemented for trees with up to 2048 leaves.
#' To analyse trees with up to 8192 leaves, you will need to a modified version
#' of the package:
#' `install.packages("BigTreeDist", repos = "https://ms609.github.io/packages/")`
#' Use `library("BigTreeDist")` *instead* of `library("TreeDist")` to load
#' the modified package &ndash; or prefix functions with the package name, e.g.
#' `BigTreeDist::TreeDistance()`.
#' 
#' As an alternative download method,
#' uninstall \pkg{TreeDist} and \pkg{TreeTools} using
#' `remove.packages()`, then use
#'  `devtools::install_github("ms609/TreeTools", ref = "more-leaves")`
#' to install the modified \pkg{TreeTools} package; then, 
#' install \pkg{TreeDist} using
#' `devtools::install_github("ms609/TreeDist", ref = "more-leaves")`.
#' (\pkg{TreeDist} will need building from source _after_ the modified 
#' \pkg{TreeTools} package has been installed, as its code links to values
#' set in the TreeTools source code.)
#' 
#' Trees with over 8192 leaves require further modification of the source code,
#' which the maintainer plans to attempt in the future; please [comment on GitHub](
#' https://github.com/ms609/TreeTools/issues/141) if you would find this useful.
#' 
#' # Parallelism
#' 
#' When `tree2 = NULL` and all trees share the same tip labels, pairwise
#' distance calculation uses a multi-threaded \acronym{OpenMP} batch path
#' automatically.  Control the number of threads with the `"mc.cores"` option:
#' 
#' ```r
#' options(mc.cores = parallel::detectCores())  # use all cores
#' options(mc.cores = 4L)                        # or a fixed number
#' ```
#' 
#' Do **not** call [`StartParallel()`] for these functions: a registered
#' R cluster disables the \acronym{OpenMP} path and replaces it with slower
#' inter-process communication.  See [`StartParallel()`] for full guidance
#' on when an R cluster is appropriate.
#' 
#' @template tree12ListParams
#' 
#' @param normalize If a numeric value is provided, this will be used as a 
#' maximum value against which to rescale results.
#' If `TRUE`, results will be rescaled against a maximum value calculated from
#' the specified tree sizes and topology, as specified in the "Normalization" 
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
#' @returns
#' If `reportMatching = FALSE`, the functions return a numeric 
#' vector specifying the requested similarities or differences.
#' 
#' If `reportMatching = TRUE`, the functions additionally return an integer
#' vector listing the index of the split in `tree2` that is matched with 
#' each split in `tree1` in the optimal matching.
#' Unmatched splits are denoted `NA`.
#' Use [`VisualizeMatching()`] to plot the optimal matching.
#' 
#' `TreeDistance()` simply returns the clustering information distance (it is
#' an alias of `ClusteringInfoDistance()`).
#'  
#' @examples 
#' tree1 <- ape::read.tree(text="((((a, b), c), d), (e, (f, (g, h))));")
#' tree2 <- ape::read.tree(text="(((a, b), (c, d)), ((e, f), (g, h)));")
#' tree3 <- ape::read.tree(text="((((h, b), c), d), (e, (f, (g, a))));")
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
#' attr(SharedPhylogeneticInfo(tree1, tree2, reportMatching = TRUE), "matching")
#' VisualizeMatching(SharedPhylogeneticInfo, tree1, tree2) # Which clades are matched?
#' 
#' DifferentPhylogeneticInfo(tree1, tree2) # Distance measure
#' DifferentPhylogeneticInfo(tree2, tree1) # The metric is symmetric
#'
#' # Are they more similar than two trees of this shape would be by chance?
#' ExpectedVariation(tree1, tree2, sample=12)["DifferentPhylogeneticInfo", "Estimate"]
#' 
#' # Every split in tree1 conflicts with every split in tree3
#' # Pairs of conflicting splits contain clustering, but not phylogenetic, 
#' # information
#' SharedPhylogeneticInfo(tree1, tree3) # = 0
#' MutualClusteringInfo(tree1, tree3) # > 0
#' 
#' # Distance functions internally convert trees to Splits objects.
#' # Pre-conversion can reduce run time if the same trees will feature in
#' # multiple comparisons
#' splits1 <- TreeTools::as.Splits(tree1)
#' splits2 <- TreeTools::as.Splits(tree2)
#' 
#' SharedPhylogeneticInfoSplits(splits1, splits2)
#' MatchingSplitInfoSplits(splits1, splits2)
#' MutualClusteringInfoSplits(splits1, splits2)
#' @template MRS 
#' 
#' @references
#' \insertAllCited{}
#' 
#' @encoding UTF-8
#' @family tree distances
#' @export
TreeDistance <- function(tree1, tree2 = NULL) {
  ClusteringInfoDistance(tree1, tree2, normalize = TRUE, reportMatching = FALSE)
}

#' @rdname TreeDistance
#' @export
SharedPhylogeneticInfo <- function(tree1, tree2 = NULL, normalize = FALSE,
                                   reportMatching = FALSE, diag = TRUE) {
  
  unnormalized <- CalculateTreeDistance(SharedPhylogeneticInfoSplits, tree1,
                                        tree2, reportMatching = reportMatching)
  
  if (diag && is.null(tree2)) {
    unnormalized <- as.matrix(unnormalized)
    diag(unnormalized) <- SplitwiseInfo(tree1)
    tree2 <- tree1
  }
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = SplitwiseInfo, Combine = .PairMean)
}

#' @rdname TreeDistance
#' @export
DifferentPhylogeneticInfo <- function(tree1, tree2 = NULL, normalize = FALSE,
                                      reportMatching = FALSE) {
  
  # Fast path (all-pairs): same tips, no matching — avoids duplicate as.Splits()
  fast <- .FastDistPath(tree1, tree2, reportMatching,
                        cpp_shared_phylo_all_pairs,
                        cpp_splitwise_info_batch)
  if (!is.null(fast)) {
    spi <- fast[["info"]]
    treesIndependentInfo <- .PairwiseSums(fast[["entropies"]])

    ret <- .FloorNumericalNoise(treesIndependentInfo - spi - spi, treesIndependentInfo)
    ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                         infoInBoth = treesIndependentInfo,
                         InfoInTree = SplitwiseInfo, Combine = "+")
    attributes(ret) <- attributes(spi)
    return(ret)
  }

  # Fast path (cross-pairs): same tips, no matching — avoids duplicate as.Splits()
  fast_many <- .FastManyManyPath(tree1, tree2, reportMatching,
                                 cpp_shared_phylo_cross_pairs,
                                 cpp_splitwise_info_batch)
  if (!is.null(fast_many)) {
    spi <- fast_many[["dists"]]
    info1 <- fast_many[["info1"]]
    info2 <- fast_many[["info2"]]
    treesIndependentInfo <- outer(info1, info2, "+")

    ret <- .FloorNumericalNoise(treesIndependentInfo - spi - spi, treesIndependentInfo)
    ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                         infoInBoth = treesIndependentInfo,
                         InfoInTree = SplitwiseInfo, Combine = "+")
    return(ret)
  }

  spi <- SharedPhylogeneticInfo(tree1, tree2, normalize = FALSE, diag = FALSE,
                                reportMatching = reportMatching)
  treesIndependentInfo <- .MaxValue(tree1, tree2, SplitwiseInfo)

  ret <- .FloorNumericalNoise(treesIndependentInfo - spi - spi, treesIndependentInfo)
  ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                       infoInBoth = treesIndependentInfo,
                       InfoInTree = SplitwiseInfo, Combine = "+")
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
ClusteringInfoDistance <- function(tree1, tree2 = NULL, normalize = FALSE,
                                   reportMatching = FALSE) {
  
  # Fast path (all-pairs): same tips, no matching — avoids duplicate as.Splits()
  fast <- .FastDistPath(tree1, tree2, reportMatching,
                        cpp_mutual_clustering_all_pairs,
                        cpp_clustering_entropy_batch)
  if (!is.null(fast)) {
    mci <- fast[["info"]]
    treesIndependentInfo <- .PairwiseSums(fast[["entropies"]])

    ret <- .FloorNumericalNoise(treesIndependentInfo - mci - mci, treesIndependentInfo)
    ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                         infoInBoth = treesIndependentInfo,
                         InfoInTree = ClusteringEntropy, Combine = "+")
    attributes(ret) <- attributes(mci)
    return(ret)
  }

  # Fast path (cross-pairs): same tips, no matching — avoids duplicate as.Splits()
  fast_many <- .FastManyManyPath(tree1, tree2, reportMatching,
                                 cpp_mutual_clustering_cross_pairs,
                                 cpp_clustering_entropy_batch)
  if (!is.null(fast_many)) {
    mci <- fast_many[["dists"]]
    info1 <- fast_many[["info1"]]
    info2 <- fast_many[["info2"]]
    treesIndependentInfo <- outer(info1, info2, "+")

    ret <- .FloorNumericalNoise(treesIndependentInfo - mci - mci, treesIndependentInfo)
    ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                         infoInBoth = treesIndependentInfo,
                         InfoInTree = ClusteringEntropy, Combine = "+")
    return(ret)
  }

  mci <- MutualClusteringInfo(tree1, tree2, normalize = FALSE, diag = FALSE,
                              reportMatching = reportMatching)
  treesIndependentInfo <- .MaxValue(tree1, tree2, ClusteringEntropy)

  ret <- .FloorNumericalNoise(treesIndependentInfo - mci - mci, treesIndependentInfo)
  ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                       infoInBoth = treesIndependentInfo,
                       InfoInTree = ClusteringEntropy, Combine = "+")
  attributes(ret) <- attributes(mci)
  
  # Return:
  ret
}

#' @export
ClusteringInfoDist <- ClusteringInfoDistance

#' @rdname TreeDistance
#' @param samples Integer specifying how many samplings to obtain; 
#' accuracy of estimate increases with `sqrt(samples)`.
#' @importFrom stats sd
#' @importFrom TreeTools as.Splits
#' @export
ExpectedVariation <- function(tree1, tree2, samples = 1e+4) {
  info1 <- SplitwiseInfo(tree1)
  info2 <- SplitwiseInfo(tree2)
  splits1 <- as.Splits(tree1)
  tipLabels <- attr(splits1, "tip.label")
  nTip <- attr(splits1, "nTip")
  splits2 <- as.Splits(tree2, tipLabels)
  
  mutualEstimates <- vapply(seq_len(samples), function(x) {
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
        "Std. Err." = ret[, "sd"] / sqrt(ret[, "n"]), 
        ret[, 2:3])
}

#' @rdname TreeDistance
#' @aliases MutualClusteringInformation
#' @export
MutualClusteringInfo <- function(tree1, tree2 = NULL, normalize = FALSE,
                                 reportMatching = FALSE, diag = TRUE) {
  unnormalized <- CalculateTreeDistance(Func = MutualClusteringInfoSplits,
                                        tree1, tree2, reportMatching)
  if (diag && is.null(tree2)) {
    unnormalized <- as.matrix(unnormalized)
    diag(unnormalized) <- ClusteringEntropy(tree1)
    tree2 <- tree1
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
SharedPhylogeneticInfoSplits <- function(splits1, splits2,
                                         nTip = attr(splits1, "nTip"),
                                         reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, nTip, cpp_shared_phylo,
                maximize = TRUE, reportMatching = reportMatching)
}

#' @rdname TreeDistance
#' @export
MutualClusteringInfoSplits <- function(splits1, splits2,
                                       nTip = attr(splits1, "nTip"),
                                       reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, nTip, cpp_mutual_clustering,
                maximize = TRUE, reportMatching = reportMatching)
}

#' Mean of two numbers
#'
#' Used for normalization and range calculation
#'
#' @keywords internal
.PairMean <- function(x, y) (x + y) / 2L

#' Largest clustering information distance reachable by one NNI move
#'
#' `NNIMaxStep()` returns the largest
#' [Clustering Information Distance][ClusteringInfoDistance] that can separate an
#' _n_-leaf tree from any tree that differs from it by a single nearest neighbour
#' interchange (NNI) move.
#'
#' A single NNI move around an internal edge exchanges two of the four subtrees
#' that meet at that edge, changing exactly one split.  Because
#' `ClusteringInfoDistance()` \insertCite{SmithDist}{TreeDist} scores trees by an
#' optimal matching of their splits, every _unchanged_ split matches its
#' counterpart perfectly and contributes nothing to the distance.  The distance of
#' an NNI move therefore depends only on the old split _S_ and the new split _S'_,
#' and hence only on the leaf counts \eqn{a, b, c, d} of the four subtrees
#' \eqn{(a + b + c + d = n)}.  Writing \eqn{H} for entropy in bits, it equals
#' their entropy distance (variation of
#' [clustering information][ClusteringEntropy], \insertCite{Meila2007}{TreeDist}):
#' \deqn{d(a, b, c, d) = 2 H(a, b, c, d) - H(a + b) - H(a + c).}
#' The maximum can thus be found deterministically, by optimizing over the four
#' subtree sizes, rather than by sampling random trees.  At the optimum the four
#' subtrees are as equal as possible (each \eqn{\lfloor n/4 \rfloor} or
#' \eqn{\lceil n/4 \rceil}), so the maximizing topology is a _local_ property of a
#' single edge, realized by _any_ tree that contains such an edge &ndash; not by a
#' globally pectinate or balanced tree.  The value is exactly two bits when four
#' divides _n_ (four equal subtrees, whose splits are then balanced and
#' independent), and slightly less otherwise.  The exchanged splits themselves
#' need not be balanced: at \eqn{n = 6} the optimal subtrees are
#' \eqn{\{2, 2, 1, 1\}} but the move runs between a \eqn{4 | 2} split and a
#' \eqn{3 | 3} split.
#'
#' @section The normalized maximum is also local:
#' With `normalize = TRUE` the distance is divided by the combined clustering
#' entropy of the two trees, \eqn{CE(T) + CE(T')}, which _does_ depend on the wider
#' topology.  The largest normalized move can nevertheless be found without
#' searching tree space, as follows.
#'
#' \eqn{T} and \eqn{T'} share every split except the one the move changes, so
#' \deqn{CE(T) + CE(T') = 2 E + H(S) + H(S'),}
#' where \eqn{E} is the summed entropy of the \eqn{n - 4} shared splits and
#' \eqn{H(S), H(S')} &ndash; the entropies of the exchanged splits &ndash; are
#' fixed by \eqn{a, b, c, d}.  Each shared split is cut by an edge lying _inside_
#' one of the four subtrees, so \eqn{E} decomposes into four independent terms,
#' \deqn{E = c(a) + c(b) + c(c) + c(d),}
#' where \eqn{c(s)} is the clustering entropy contributed by a subtree of \eqn{s}
#' leaves.  The numerator \eqn{d(a, b, c, d)} is independent of the subtree
#' _shapes_, so for fixed sizes the ratio is largest when each \eqn{c(s)} is as
#' _small_ as possible &ndash; and each subtree is minimized independently.  The
#' minimum \eqn{c(s)} satisfies the recursion
#' \deqn{c(1) = 0, \qquad c(s) = H(s) + \min_{1 \le i \le s/2}\{c(i) + c(s - i)\},}
#' computed once by dynamic programming.  The largest normalized NNI move is then
#' \deqn{\max_{a, b, c, d} \frac{d(a, b, c, d)}{2[c(a) + c(b) + c(c) + c(d)] +
#' H(S) + H(S')},}
#' a purely local optimization over the four subtree sizes.  (The minimizing
#' subtree shape is not in general balanced; for \eqn{s = 6}, for instance, it
#' splits \eqn{2 + 4} rather than \eqn{3 + 3}.)
#'
#' @param tree Object of supported class representing a tree or list of trees,
#' or an integer specifying the number of leaves in a tree/trees.
#' @param normalize Logical specifying whether to normalize the distance against
#' the summed clustering information of the two trees, giving the largest
#' _normalized_ Clustering Information Distance attainable by one NNI move.
#'
#' @return `NNIMaxStep()` returns a numeric vector, one entry per tree (or leaf
#' count), giving the largest attainable distance &ndash; in bits when
#' `normalize = FALSE`, or as a fraction in the range \[0, 1] when
#' `normalize = TRUE`.  Entries are `NA` where _n_ &lt; 4, as no NNI move exists.
#' Two attributes record the maximizing configuration: `"subtrees"`, the sizes of
#' the four subtrees around the moved edge, and `"splits"`, the sizes of the two
#' splits that the move exchanges.  Where more than one tree is supplied, each
#' attribute is a list with one entry per tree.
#'
#' @examples
#' # Largest clustering information distance from a single NNI move
#' NNIMaxStep(8)  # exactly two bits: eight is a multiple of four
#' NNIMaxStep(6)  # a little less
#'
#' # Read off the maximizing local topology
#' m6 <- NNIMaxStep(6)
#' attr(m6, "subtrees")
#' attr(m6, "splits")
#'
#' # Vectorized over leaf counts, and accepting a tree
#' NNIMaxStep(4:12)
#' library("TreeTools", quietly = TRUE)
#' NNIMaxStep(BalancedTree(19))
#'
#' # Normalized: solved locally, without searching tree space
#' NNIMaxStep(12, normalize = TRUE)
#'
#' @template MRS
#' @family tree distances
#' @seealso
#' The distance itself: [`ClusteringInfoDistance()`]
#'
#' Diameter of the NNI metric: [`NNIDiameter()`]
#' @references \insertAllCited{}
#' @export
NNIMaxStep <- function(tree, normalize = FALSE) {
  UseMethod("NNIMaxStep")
}

#' @export
NNIMaxStep.numeric <- function(tree, normalize = FALSE) {
  res <- lapply(tree, .NNIMaxStep1, normalize = normalize)
  value <- vapply(res, `[[`, double(1), "value")
  subtrees <- lapply(res, `[[`, "subtrees")
  splits <- lapply(res, `[[`, "splits")
  if (length(tree) == 1L) {
    attr(value, "subtrees") <- subtrees[[1]]
    attr(value, "splits") <- splits[[1]]
  } else {
    attr(value, "subtrees") <- subtrees
    attr(value, "splits") <- splits
  }
  value
}

#' @importFrom TreeTools NTip
#' @export
NNIMaxStep.phylo <- function(tree, normalize = FALSE) {
  NNIMaxStep(NTip(tree), normalize = normalize)
}

#' @export
NNIMaxStep.multiPhylo <- NNIMaxStep.phylo

#' @export
NNIMaxStep.list <- function(tree, normalize = FALSE) {
  lapply(tree, NNIMaxStep, normalize = normalize)
}

# Binary split entropy of a k | (n - k) split, in bits (vectorized over k).
.BinaryEntropyBits <- function(k, n) {
  p <- k / n
  ifelse(p <= 0 | p >= 1, 0, -(p * log2(p) + (1 - p) * log2(1 - p)))
}

# Minimum clustering-entropy contribution (bits) of a rooted subtree of each
# size 1..n, splits measured against the whole tree's n leaves. A subtree of
# size `s` contributes its attaching split H(s, n) plus, recursively, its two
# daughters; the daughter split is chosen to minimize the total. `split[s]`
# records the smaller daughter's size, so the subtree shape is recoverable.
.MinSubtreeEntropy <- function(n) {
  cost <- numeric(n)
  split <- integer(n)
  if (n >= 2L) {
    cost[2] <- .BinaryEntropyBits(2L, n)
    split[2] <- 1L
  }
  for (s in seq_len(n)[-(1:2)]) {
    i <- seq_len(s %/% 2L)
    daughters <- cost[i] + cost[s - i]
    best <- which.min(daughters)
    cost[s] <- .BinaryEntropyBits(s, n) + daughters[best]
    split[s] <- i[best]
  }
  list(cost = cost, split = split)
}

# Largest NNI clustering information distance for a single leaf count.
.NNIMaxStep1 <- function(n, normalize) {
  n <- as.integer(round(n))
  if (is.na(n) || n < 4L) {
    return(list(value = NA_real_, subtrees = NULL, splits = NULL))
  }
  cost <- if (normalize) .MinSubtreeEntropy(n)[["cost"]] else NULL
  gLog <- function(x) (x / n) * log2(x / n)  # -contribution to joint entropy

  best <- -Inf
  bestParts <- NULL
  bestSplits <- NULL
  # Enumerate every partition of n into four parts p <= q <= r <= s: each is the
  # multiset of subtree sizes around some edge, and every such multiset is
  # realizable, so this is exhaustive over all (tree, NNI move) pairs.
  for (p in seq_len(n %/% 4L)) {
    for (q in p:((n - p) %/% 3L)) {
      # q <= (n - p) %/% 3 guarantees rMax >= q, so r <- q:rMax is non-empty.
      r <- q:((n - p - q) %/% 2L)
      s <- n - p - q - r
      # The four subtrees define three possible splits, with one side p + q,
      # p + r, p + s. An NNI move exchanges two of them; the largest distance
      # keeps the two least balanced (smallest entropy), i.e. drops the largest.
      h1 <- .BinaryEntropyBits(p + q, n)
      h2 <- .BinaryEntropyBits(p + r, n)
      h3 <- .BinaryEntropyBits(p + s, n)
      keptEntropy <- h1 + h2 + h3 - pmax(h1, h2, h3)
      joint <- -(gLog(p) + gLog(q) + gLog(r) + gLog(s))
      vi <- 2 * joint - keptEntropy
      value <- if (normalize) {
        vi / (2 * (cost[p] + cost[q] + cost[r] + cost[s]) + keptEntropy)
      } else {
        vi
      }
      winner <- which.max(value)
      if (value[winner] > best) {
        best <- value[winner]
        rw <- r[winner]
        sw <- s[winner]
        bestParts <- c(p, q, rw, sw)
        sizes <- c(p + q, p + rw, p + sw)
        entropies <- c(h1, h2[winner], h3[winner])
        bestSplits <- sort(sizes[-which.max(entropies)])
      }
    }
  }
  list(value = best, subtrees = bestParts, splits = bestSplits)
}
