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
PartitionInfo <- function(tree) {
  if (class(tree) == 'phylo') {
    PartitionInfoSplits(Tree2Splits(tree))
  } else {
    splits <- lapply(tree, Tree2Splits)
    vapply(splits, PartitionInfoSplits, double(1L))
  }
}

#' @describeIn PartitionInfo Takes splits instead of trees
#' @export
PartitionInfoSplits <- function(splits) {
  nTerminals <- nrow(splits)
  inSplit <- colSums(splits)
  
  sum(vapply(inSplit, LnRooted.int, 0) + 
        + vapply(nTerminals - inSplit,  LnRooted.int, 0)
      - LnUnrooted.int(nTerminals)
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
#' @references \insertRef{Meila2007}{TreeSearch}
#' @author Martin R. Smith
#' @keywords internal
#' @export
ClusteringInfo <- function(tree) {
  if (class(tree) == 'phylo') {
    ClusteringInfoSplits(Tree2Splits(tree))
  } else {
    splits <- lapply(tree, Tree2Splits)
    vapply(splits, ClusteringInfoSplits, double(1L))
  }
}

#' @describeIn ClusteringInfo Takes splits instead of trees
#' @export
ClusteringInfoSplits <- function (splits) {
  nTip <- nrow(splits)
  inSplit <- colSums(splits)
  splitP <- rbind(inSplit, nTip - inSplit) / nTip
  
  # Entropy measures the bits required to transmit the cluster label of each tip.
  # The total information content is thus entropy * nTip.
  # See Vinh (2010: p. 2840) 
  # 
  # Return:
  sum(apply(splitP, 2, Entropy)) / log(2) * nTip
}
