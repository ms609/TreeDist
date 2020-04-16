#' Information content of splits within a tree
#' 
#' Sums the phylogenetic information content for all splits within a 
#' phylogenetic tree.  This value will be greater than the total information 
#' content of the tree where a tree contains multiple splits, as 
#' these splits will contain mutual information.
#' 
#' @param x A tree of class `phylo`, a list of trees, or a `multiPhylo` object.
#' 
#' @family information functions
#' 
#' @seealso 
#' 
#' An introduction to the phylogenetic information content of a split is given
#' in [`SplitInformation()`](https://ms609.github.io/TreeTools/reference/SplitInformation.html)
#' and in a [package vignette](https://ms609.github.io/TreeDist/articles/information.html).
#' 
#' @template MRS
#' @export
SplitwiseInfo <- function (x) UseMethod('SplitwiseInfo')

#' @export
SplitwiseInfo.phylo <- function (x) SplitwiseInfo.Splits(as.Splits(x))
  
#' @export
SplitwiseInfo.multiPhylo <- function (x) {
  vapply(as.Splits(x), SplitwiseInfo, double(1L))
}

#' @export
SplitwiseInfo.list <- SplitwiseInfo.multiPhylo

#' @importFrom TreeTools LnRooted.int LnUnrooted.int TipsInSplits
#' @export
SplitwiseInfo.Splits <- function(x) {
  nTip <- attr(x, 'nTip')
  inSplit <- TipsInSplits(x)
  
  sum(vapply(inSplit, LnRooted.int, 0) + 
        + vapply(nTip - inSplit,  LnRooted.int, 0)
      - LnUnrooted.int(nTip)
  ) / -log(2)
}

#' Clustering entropy of all splits within a tree
#' 
#' Sums the entropy (`ClusteringEntropy()`) or information content
#' (`ClusteringInfo()`) across each split within a phylogenetic tree, treating
#' each split as dividing the leaves of the tree into two clusters (_sensu_
#' Meila 2007; Vinh _et al._ 2010).
#' 
#' As entropy measures the bits required to transmit the cluster label of each 
#' leaf (Vinh _et al._ 2010: p. 2840), the information content of a split is its
#' entropy multiplied by the number of leaves. 
#' 
#' @param x A tree of class `phylo`, a list of trees, or a `multiPhylo` object.
#' 
#' @return The sum of the entropies or (clustering) information content, in bits,
#' of each split in `x`.
#' 
#' @examples 
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
#' @references
#' 
#' - \insertRef{Meila2007}{TreeDist}
#' 
#' - \insertRef{Vinh2010}{TreeDist}
#' 
#' 
#' @family information functions
#' @template MRS
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
