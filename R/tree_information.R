#' Information content of partitions within a tree
#' 
#' Sums the phylogenetic information content for all partitions within a 
#' phylogenetic tree.  This value will be greater than the total information 
#' content of the tree where a tree contains multiple partitions, as 
#' these partitions will contain mutual information
#' 
#' @param tree A tree of class `phylo`, a list of trees, or a `multiPhylo` object.
#' 
#' @author Martin R. Smith
#' @keywords internal
#' @export
PartitionInfo <- function (x) UseMethod('PartitionInfo')

#' @export
PartitionInfo.phylo <- function (x) PartitionInfo.Splits(as.Splits(x))
  
#' @export
PartitionInfo.list <- PartitionInfo.multiPhylo <- function (x) {
  vapply(as.Splits(x), PartitionInfo, double(1L))
}

#' @importFrom TreeTools LnRooted.int LnUnrooted.int TipsInSplits
#' @export
PartitionInfo.Splits <- function(x) {
  nTip <- attr(x, 'nTip')
  inSplit <- TipsInSplits(x)
  
  sum(vapply(inSplit, LnRooted.int, 0) + 
        + vapply(nTip - inSplit,  LnRooted.int, 0)
      - LnUnrooted.int(nTip)
  ) / -log(2)
}


#' Clustering information content of all partitions within a tree
#' 
#' Sums the clustering information content (Meila 2007) across each partition within a 
#' phylogenetic tree.  This value will be greater than the total information 
#' content of the tree where a tree contains multiple partitions, as 
#' these partitions will contain mutual information.
#' 
#' Note that Meila (2007) and Vinh _et al_. (2010) deal with entropy, which 
#' denotes the bits necessary to denote the cluster to which each tip belongs,
#' i.e. bits/tip (see Vinh _et al._ 2010).
#' This function multiplies the entropy by the number of tips to produce a
#' measure of the information content of the tree as a whole.
#' 
#' @param tree A tree of class `phylo`, a list of trees, or a `multiPhylo` object.
#' 
#' @return Information content, in bits.
#' 
#' @references \insertRef{Meila2007}{TreeDist}
#' @author Martin R. Smith
#' @keywords internal
#' @importFrom TreeTools as.Splits
#' @export
ClusteringInfo <- function (x) UseMethod("ClusteringInfo")

#' @export
ClusteringInfo.phylo <- function (x) ClusteringInfo.Splits(as.Splits(x))

#' @export
ClusteringInfo.list <- ClusteringInfo.multiPhylo <- function (x)
    vapply(as.Splits(x), ClusteringInfo.Splits, double(1L))

#' @export
ClusteringInfo.Splits <- function (x) {
  nTip <- attr(x, 'nTip')
  inSplit <- TipsInSplits(x)
  splitP <- rbind(inSplit, nTip - inSplit) / nTip
  
  # Entropy measures the bits required to transmit the cluster label of each tip.
  # The total information content is thus entropy * nTip.
  # See Vinh (2010: p. 2840) 
  # 
  # Return:
  sum(apply(splitP, 2, Entropy)) * nTip
}
