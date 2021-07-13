#' Information content of splits within a tree
#' 
#' Sum the entropy (`ClusteringEntropy()`), clustering information content
#' (`ClusteringInfo()`), or phylogenetic information content (`SplitwiseInfo()`)
#' across each split within a phylogenetic tree,
#' or the consensus of a set of phylogenetic trees (`ConsensusInfo()`).
#' This value will be greater than the total information 
#' content of the tree where a tree contains multiple splits, as 
#' these splits are not independent and thus contain mutual information that is
#' counted twice.
#' 
#' @section Clustering information:
#' 
#' Clustering entropy addresses the question "how much information is contained
#' in the splits within a tree". Its approach is complementary to the 
#' phylogenetic information content, used in [`SplitwiseInfo()`].
#' In essence, it asks, given a split that subdivides the leaves of a tree into
#' two partitions, how easy it is to predict which partition a randomly drawn 
#' leaf belongs to \insertCite{@Meila2007; @Vinh2010}{TreeDist}.
#' 
#' 
#' Formally, the entropy of a split _S_ that divides _n_ leaves into two
#' partitions of sizes _a_ and _b_ is given by 
#' _H(S)_ = - _a/n_ log _a/n_ - _b/n_ log _b/n_.
#' 
#' Base 2 logarithms are conventionally used, such that entropy is measured in
#' bits.  
#' Entropy denotes the number of bits that are necessary to encode the outcome
#' of a random variable: here, the random variable is "what partition does a
#' randomly selected leaf belong to". 
#' 
#' An even split has an entropy of 1 bit: there is no better way of encoding 
#' an outcome than using one bit to specify which of the two partitions the 
#' randomly selected leaf belongs to.
#' 
#' An uneven split has a lower entropy: membership of the larger partition is
#' common, and thus less surprising; it can be signified using fewer bits in an
#' optimal compression system.
#' 
#' If this sounds confusing, let's consider creating a code to transmit the 
#' cluster label of two randomly selected leaves.  One straightforward
#' option would be to use 
#' 
#' - `00` = 'Both leaves belong to partition A'
#' 
#' - `11` = 'Both leaves belong to partition B'
#' 
#' - `01` = 'First leaf in A, second in B`
#' 
#' - `10` = 'First leaf in B, second in A`
#' 
#' This code uses two bits to transmit the partition labels of two leaves.
#' If partitions A and B are equiprobable, this is the optimal code; our 
#' entropy -- the average information content required per leaf -- is 1 bit.
#' 
#' 
#' Alternatively, we could use the (suboptimal) code
#' 
#' - `0` = 'Both leaves belong to partition A'
#' 
#' - `111` = 'Both leaves belong to partition B'
#' 
#' - `101` = 'First leaf in A, second in B`
#' 
#' - `110` = 'First leaf in B, second in A`
#' 
#' If A is much larger than B, then most pairs of leaves will require just
#' a single bit (code `0`). The additional bits when 1+ leaf belongs to B
#' may be required sufficiently rarely that the average message 
#' requires fewer than two bits for two leaves, so the entropy is less than 
#' 1 bit.  (The optimal coding strategy will depend on the exact sizes
#' of A and B.)
#' 
#' 
#' As entropy measures the bits required to transmit the cluster label of each 
#' leaf \insertCite{@@Vinh2010: p. 2840}{TreeDist}, the information content of 
#' a split is its entropy multiplied by the number of leaves. 
#' 
#' @section Phylogenetic information:
#' 
#' Phylogenetic information expresses the information content of a split
#' in terms of the probability that a uniformly selected tree will contain it
#' \insertCite{Thorley1998}{TreeDist}.
#' 
#' @section Consensus information:
#' 
#' The information content of splits in a consensus tree is calculated by
#' interpreting support values (i.e. the proportion of trees containing
#' each split in the consensus) as probabilities that the true tree contains
#' that split, following \insertCite{SmithCons;textual}{TreeDist}.
#' 
#' @return `SplitwiseInfo()`, `ClusteringInfo()` and `ClusteringEntropy()`
#' return the splitwise information content of the tree -- or of each split
#' in turn, if `sum = FALSE` -- in bits.
#' 
#' @return `ConsensusInfo()` returns the splitwise information content of the
#' majority rule consensus of `trees`.
#' 
#' @seealso 
#' 
#' An introduction to the phylogenetic information content of a split is given
#' in [`SplitInformation()`](https://ms609.github.io/TreeTools/reference/SplitInformation.html)
#' and in a [package vignette](https://ms609.github.io/TreeDist/articles/information.html).
#' 
#' @examples
#' library("TreeTools", quietly = TRUE)
#' 
#' SplitwiseInfo(PectinateTree(8))
#' tree <- read.tree(text = "(a, b, (c, (d, e, (f, g)0.8))0.9);")
#' SplitwiseInfo(tree)
#' SplitwiseInfo(tree, TRUE)
#' @template MRS
#' 
#' @references
#' \insertAllCited{}
#' 
#' @family information functions
#' @encoding UTF-8
#' @name TreeInfo
NULL


#' @param x A tree of class `phylo`, a list of trees, or a `multiPhylo` object.
#' @param p A vector of probabilities corresponding to each split in `x`. 
#' Specify `TRUE` to calculate this vector using the node labels of each tree.
#' @param sum Logical: if `TRUE`, sum the information content of each split to
#' provide the total splitwise information content of the tree.
#' @rdname TreeInfo
#' @export
SplitwiseInfo <- function (x, p = NULL, sum = TRUE) {
  UseMethod('SplitwiseInfo')
}

.GetPFromLabels <- function (tree, p, splits = as.Splits(tree)) {
  if (length(p) == 1L) { # length(NULL) == 0
    if (p == FALSE) {
      p <- NULL
    } else {
      np <- tree$node.label[as.integer(names(splits)) - NTip(tree)]
      if (is.null(np)) {
        np <- rep_len(p, length(splits))
      }
      p <- as.double(np) / p
      p[is.na(p)] <- 1
      if (any(p > 1)) {
        stop("Nodes must be labelled with probabilities <= 1")
      }
    }
  }
  # Return:
  p
}

#' @export
SplitwiseInfo.phylo <- function (x, p = NULL, sum = TRUE) {
  splits <- as.Splits(x)
  SplitwiseInfo.Splits(splits, .GetPFromLabels(x, p, splits), sum)
}
  
#' @export
SplitwiseInfo.multiPhylo <- function (x, p = NULL, sum = TRUE) {
  vapply(as.Splits(x), SplitwiseInfo, double(1L), p, sum)
}

#' @export
SplitwiseInfo.list <- SplitwiseInfo.multiPhylo

#' @importFrom TreeTools Log2Rooted.int Log2Unrooted.int TipsInSplits
#' @export
SplitwiseInfo.Splits <- function(x, p = NULL, sum = TRUE) {
  nTip <- attr(x, 'nTip')
  inSplit <- TipsInSplits(x)
  
  if (is.null(p)) {
    ret <- Log2Unrooted.int(nTip) -
      vapply(inSplit, Log2Rooted.int, 0) -
      vapply(nTip - inSplit, Log2Rooted.int, 0)
    if (sum) sum(ret) else ret
  } else {
    q <- 1L - p
    qNonZero <- as.logical(1L - p)
    q <- q[qNonZero]
    
    l2n <- Log2Unrooted.int(nTip)
    
    l2nConsistent <- vapply(inSplit, Log2Rooted.int, 0) +
      vapply(nTip - inSplit, Log2Rooted.int, 0)
    
    l2pConsistent <- l2nConsistent - l2n
    l2pInconsistent <- log2(-expm1(l2pConsistent[qNonZero] * log(2)))
    
    l2nInconsistent <- l2pInconsistent + l2n
    
    # Return:
    if (sum) {
      sum(l2n * length(inSplit),
          p * (log2(p) - l2nConsistent),
          q * (log2(q) - l2nInconsistent))
    } else {
      ret <- l2n + (p * (log2(p) - l2nConsistent))
      ret[names(l2nInconsistent)] <- ret[names(l2nInconsistent)] +
        (q * (log2(q) - l2nInconsistent))
      ret
    }
  }
}

#' @export
SplitwiseInfo.NULL <- function (x, p = NULL, sum = TRUE) NULL

#' @examples
#' 
#' # Clustering entropy of an even split = 1 bit
#' ClusteringEntropy(TreeTools::as.Splits(c(rep(TRUE, 4), rep(FALSE, 4))))
#' 
#' # Clustering entropy of an uneven split < 1 bit
#' ClusteringEntropy(TreeTools::as.Splits(c(rep(TRUE, 2), rep(FALSE, 6))))
#' 
#' tree1 <- TreeTools::BalancedTree(8)
#' tree2 <- TreeTools::PectinateTree(8)
#' 
#' ClusteringInfo(tree1)
#' ClusteringEntropy(tree1)
#' ClusteringInfo(list(one = tree1, two = tree2))
#' 
#' ClusteringInfo(tree1) + ClusteringInfo(tree2)
#' ClusteringEntropy(tree1) + ClusteringEntropy(tree2)
#' ClusteringInfoDistance(tree1, tree2)
#' MutualClusteringInfo(tree1, tree2)
#' 
#' # Clustering entropy with uncertain splits
#' tree <- ape::read.tree(text = "(a, b, (c, (d, e, (f, g)0.8))0.9);")
#' ClusteringInfo(tree)
#' ClusteringInfo(tree, TRUE)
#' @importFrom TreeTools as.Splits
#' @rdname TreeInfo
#' @export
ClusteringEntropy <- function (x, p = NULL, sum = TRUE) {
  UseMethod("ClusteringEntropy")
}

#' @rdname TreeInfo
#' @export
ClusteringInfo <- function (x, p = NULL, sum = TRUE) {
  UseMethod("ClusteringInfo")
}

#' @rdname TreeInfo
#' @export
ClusteringEntropy.phylo <- function (x, p = NULL, sum = TRUE) {
  splits <- as.Splits(x)
  ClusteringEntropy.Splits(splits, p = .GetPFromLabels(x, p, splits), sum)
}

#' @rdname TreeInfo
#' @export
ClusteringEntropy.list <- function (x, p = NULL, sum = TRUE)
    vapply(as.Splits(x), ClusteringEntropy.Splits, double(1L), p = p, sum)

#' @rdname TreeInfo
#' @export
ClusteringEntropy.multiPhylo <- ClusteringEntropy.list

#' @rdname TreeInfo
#' @export
ClusteringEntropy.Splits <- function (x, p = NULL, sum = TRUE) {
  nLeaves <- attr(x, 'nTip')
  inSplit <- TipsInSplits(x)
  splitP <- rbind(inSplit, nLeaves - inSplit, deparse.level = 0L) / nLeaves
  if (is.null(p)) {
    p <- 1L
  }
  
  ret <- p * apply(splitP, 2, Entropy)
  # Return:
  if (sum) sum(ret) else ret
}

#' @export
ClusteringEntropy.NULL <- function (x, p = NULL, sum = TRUE) NULL

#' @rdname TreeInfo
#' @export
ClusteringInfo.phylo <- function (x, p = NULL, sum = TRUE) {
  splits <- as.Splits(x)
  ClusteringInfo.Splits(splits, p = .GetPFromLabels(x, p, splits), sum)
}

#' @rdname TreeInfo
#' @export
ClusteringInfo.list <- function (x, p = NULL, sum = TRUE)
    vapply(as.Splits(x), ClusteringInfo.Splits, double(1L), p = p, sum)

#' @rdname TreeInfo
#' @export
ClusteringInfo.multiPhylo <- ClusteringInfo.list

#' @rdname TreeInfo
#' @export
ClusteringInfo.Splits <- function (x, p = NULL, sum = TRUE) {
  nLeaves <- attr(x, 'nTip')
  ClusteringEntropy.Splits(x, p, sum) * nLeaves
}

#' @export
ClusteringInfo.NULL <- function (x, p = NULL, sum = TRUE) NULL

#' @rdname TreeInfo
#' 
#' @param trees List of `phylo` objects, optionally with class `multiPhylo`.
#' @param info Abbreviation of 'phylogenetic' or 'clustering', specifying
#' the concept of information to employ.
#' @param check.tips Logical specifying whether to renumber leaves such that
#' leaf numbering is consistent in all trees.
#' 
#' @examples
#' 
#' # Support-weighted information content of a consensus tree
#' set.seed(0)
#' trees <- list(RandomTree(8), RootTree(BalancedTree(8), 1), PectinateTree(8))
#' cons <- consensus(trees, p = 0.5)
#' p <- SplitFrequency(cons, trees) / length(trees)
#' plot(cons)
#' LabelSplits(cons, signif(SplitwiseInfo(cons, p, sum = FALSE), 4))
#' ConsensusInfo(trees)
#' LabelSplits(cons, signif(ClusteringInfo(cons, p, sum = FALSE), 4))
#' ConsensusInfo(trees, 'clustering')
#' @export
ConsensusInfo <- function (trees, info = 'phylogenetic', check.tips = TRUE) {
  mode <- pmatch(tolower(info),
                 c('phylogenetic', 'clustering', 'spic', 'scic')) %% 2
  if (is.na(mode)) {
    stop("`info` must be 'phylogenetic' or 'clustering'")
  }
  if (inherits(trees, 'phylo')) {
    # Convert to length-1 multiphylo object
    trees <- c(trees)
  }
  if (length(trees) == 1L) {
    return((if (mode) SplitwiseInfo else ClusteringInfo)(trees))
  }
  if (check.tips) {
    trees <- RenumberTips(trees, trees[[1]])
  }
  consensus_info(trees, mode == 1L)
}
