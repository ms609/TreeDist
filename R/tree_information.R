#' Information content of splits within a tree
#' 
#' Sum the phylogenetic information content for all splits within a 
#' phylogenetic tree.  This value will be greater than the total information 
#' content of the tree where a tree contains multiple splits, as 
#' these splits will contain mutual information.
#' 
#' @param x A tree of class `phylo`, a list of trees, or a `multiPhylo` object.
#' @param p A vector of probabilities corresponding t oeach split in `x`. 
#' Specify `TRUE` to calculate this vector using the node labels of each tree.
#' 
#' @family information functions
#' 
#' @seealso 
#' 
#' An introduction to the phylogenetic information content of a split is given
#' in [`SplitInformation()`](https://ms609.github.io/TreeTools/reference/SplitInformation.html)
#' and in a [package vignette](https://ms609.github.io/TreeDist/articles/information.html).
#' 
#' @examples
#' SplitwiseInfo(TreeTools::PectinateTree(8))
#' tree <- ape::read.tree(text = "(a, b, (c, (d, e, (f, g)0.8))0.9);")
#' SplitwiseInfo(tree)
#' SplitwiseInfo(tree, TRUE)
#' @template MRS
#' @export
SplitwiseInfo <- function (x, p = NULL) UseMethod('SplitwiseInfo')

#' @export
SplitwiseInfo.phylo <- function (x, p = NULL) {
  splits <- as.Splits(x)
  if (length(p) == 1L) {
    np <- x$node.label[as.integer(names(splits)) - NTip(x)]
    p <- as.double(np) / p
    p[is.na(p)] <- 1
    if (any(p > 1)) {
      stop("Nodes must be labelled with probabilities <= 1")
    }
  }
  SplitwiseInfo.Splits(splits, p)
}
  
#' @export
SplitwiseInfo.multiPhylo <- function (x, p = NULL) {
  vapply(as.Splits(x), SplitwiseInfo, double(1L), p)
}

#' @export
SplitwiseInfo.list <- SplitwiseInfo.multiPhylo

#' @importFrom TreeTools Log2Rooted.int Log2Unrooted.int TipsInSplits
#' @export
SplitwiseInfo.Splits <- function(x, p = NULL) {
  nTip <- attr(x, 'nTip')
  inSplit <- TipsInSplits(x)
  
  if (is.null(p)) {
    sum(Log2Unrooted.int(nTip) -
           vapply(inSplit, Log2Rooted.int, 0) -
           vapply(nTip - inSplit, Log2Rooted.int, 0)
    )
  } else {
    
    #p <- 0.6 * c(1, 0, 0) + 0.4 * c(0, 1/2, 1/2)
    
    # expect 1.37
    
    q <- 1L - p
    qNonZero <- as.logical(1L - p)
    
    originalH <- Log2Unrooted.int(nTip)
    l2StartP <- -originalH
    
    l2pConsistent <- -(Log2Unrooted.int(nTip) -
      vapply(inSplit, Log2Rooted.int, 0) -
      vapply(nTip - inSplit, Log2Rooted.int, 0))
    
    l2pInconsistent <- log2(-expm1(l2pConsistent[qNonZero] * log(2)))
           
    log2p <- log2(p) + (l2StartP - l2pConsistent)
    log2q <- log2(q[qNonZero]) + (l2StartP - l2pInconsistent)
    
    # Return:
    -sum(-(p * log2p) - (q * log2q))
  }
}

#' @export
SplitwiseInfo.NULL <- function (x, p = NULL) 0

#' Clustering entropy of all splits within a tree
#' 
#' Sum the entropy (`ClusteringEntropy()`) or information content
#' (`ClusteringInfo()`) across each split within a phylogenetic tree, treating
#' each split as dividing the leaves of the tree into two clusters (_sensu_
#' Meil\ifelse{html}{\out{&#259;}}{a} 2007; Vinh _et al._ 2010).
#' 
#' Clustering entropy addresses the question "how much information is contained
#' in the splits within a tree". Its approach is complementary to the 
#' phylogenetic information content, used in [`SplitwiseInfo()`].
#' In essence, it asks, given a split that subdivides the leaves of a tree into
#' two partitions, how easy it is to predict which partition a randomly drawn 
#' leaf belongs to.
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
#' leaf (Vinh _et al._ 2010: p. 2840), the information content of a split is its
#' entropy multiplied by the number of leaves. 
#' 
#' @param x A tree of class `phylo`, a list of trees, or a `multiPhylo` object.
#' 
#' @return Returns the sum of the entropies or (clustering) information content, in bits,
#' of each split in `x`.
#' 
#' @examples
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
#' @template MRS
#' 
#' @references
#' 
#' - \insertRef{Meila2007}{TreeDist}
#' 
#' - \insertRef{Vinh2010}{TreeDist}
#' 
#' 
#' @encoding UTF-8
#' @family information functions
#' @importFrom TreeTools as.Splits
#' @export
ClusteringEntropy <- function (x) UseMethod("ClusteringEntropy")

#' @rdname ClusteringEntropy
#' @export
ClusteringInfo <- function (x) UseMethod("ClusteringInfo")

#' @rdname ClusteringEntropy
#' @export
ClusteringEntropy.phylo <- function (x) ClusteringEntropy.Splits(as.Splits(x))

#' @rdname ClusteringEntropy
#' @export
ClusteringEntropy.list <- function (x)
    vapply(as.Splits(x), ClusteringEntropy.Splits, double(1L))

#' @rdname ClusteringEntropy
#' @export
ClusteringEntropy.multiPhylo <- ClusteringEntropy.list

#' @rdname ClusteringEntropy
#' @export
ClusteringEntropy.Splits <- function (x) {
  nLeaves <- attr(x, 'nTip')
  inSplit <- TipsInSplits(x)
  splitP <- rbind(inSplit, nLeaves - inSplit, deparse.level = 0L) / nLeaves
  
  # Return:
  sum(apply(splitP, 2, Entropy))
}

#' @export
ClusteringEntropy.NULL <- function (x) NULL

#' @rdname ClusteringEntropy
#' @export
ClusteringInfo.phylo <- function (x) ClusteringInfo.Splits(as.Splits(x))

#' @rdname ClusteringEntropy
#' @export
ClusteringInfo.list <- function (x)
    vapply(as.Splits(x), ClusteringInfo.Splits, double(1L))

#' @rdname ClusteringEntropy
#' @export
ClusteringInfo.multiPhylo <- ClusteringInfo.list

#' @rdname ClusteringEntropy
#' @export
ClusteringInfo.Splits <- function (x) {
  nLeaves <- attr(x, 'nTip')
  inSplit <- TipsInSplits(x)
  splitP <- rbind(inSplit, nLeaves - inSplit, deparse.level = 0L) / nLeaves
  
  # Return:
  sum(apply(splitP, 2, Entropy)) * nLeaves
}

#' @export
ClusteringInfo.NULL <- function (x) NULL
