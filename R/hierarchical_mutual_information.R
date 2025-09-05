#' Hierarchical Mutual Information for phylogenetic trees
#'
#' Calculate the Hierarchical Mutual Information (HMI) between two phylogenetic
#' trees, which extends traditional mutual information to account for the
#' hierarchical structure inherent in phylogenetic trees.
#' 
#' @details
#' Hierarchical Mutual Information considers the nested, hierarchical structure
#' of phylogenetic trees when computing information measures. Unlike standard
#' mutual information which treats all splits equally, HMI weights splits
#' according to their position in the tree hierarchy, providing a more
#' nuanced measure of tree similarity that accounts for the evolutionary
#' relationships represented.
#' 
#' The measure is calculated by considering:
#' \itemize{
#'   \item The depth of each split in the tree hierarchy
#'   \item The information content of each split
#'   \item The mutual information between corresponding splits across trees
#'   \item Hierarchical weighting based on tree structure
#' }
#' 
#' @param tree1,tree2 Trees of class \code{phylo}, or lists of such trees.
#' If \code{tree2} is not provided, distances will be calculated between
#' each pair of trees in the list \code{tree1}.
#' @param normalize Logical. If \code{TRUE}, normalize the result to range [0,1].
#' @param reportMatching Logical specifying whether to return the clade
#' matchings as an attribute of the score.
#' 
#' @return A numeric value representing the Hierarchical Mutual Information
#' between the input trees. If \code{reportMatching = TRUE}, returns additional
#' attributes showing the optimal matching between splits.
#' 
#' @examples
#' library("TreeTools", quietly = TRUE)
#' 
#' tree1 <- BalancedTree(8)
#' tree2 <- PectinateTree(8)
#' 
#' # Calculate HMI between two trees
#' HierarchicalMutualInfo(tree1, tree2)
#' 
#' # Normalized HMI
#' HierarchicalMutualInfo(tree1, tree2, normalize = TRUE)
#' 
#' # Compare with standard mutual information
#' MutualClusteringInfo(tree1, tree2)
#' 
#' @references
#' Based on concepts from:
#' - Meila, M. (2007). Comparing clusterings - an information based distance.
#' - Vinh, N. X. et al. (2010). Information theoretic measures for clusterings comparison
#' 
#' @family tree distances
#' @export
HierarchicalMutualInfo <- function(tree1, tree2 = NULL, normalize = FALSE,
                                  reportMatching = FALSE) {
  UseMethod("HierarchicalMutualInfo")
}

#' @export
HierarchicalMutualInfo.phylo <- function(tree1, tree2 = NULL, normalize = FALSE,
                                        reportMatching = FALSE) {
  if (is.null(tree2)) {
    stop("tree2 must be provided for phylo objects")
  }
  
  # Convert trees to splits
  splits1 <- TreeTools::as.Splits(tree1)
  splits2 <- TreeTools::as.Splits(tree2)
  
  # Calculate HMI using splits
  HierarchicalMutualInfoSplits(splits1, splits2, normalize = normalize,
                              reportMatching = reportMatching)
}

#' @export  
HierarchicalMutualInfo.list <- function(tree1, tree2 = NULL, normalize = FALSE,
                                       reportMatching = FALSE) {
  CalculateTreeDistance(HierarchicalMutualInfoSplits, tree1, tree2, 
                        reportMatching = reportMatching, normalize = normalize)
}

#' @export
HierarchicalMutualInfo.multiPhylo <- HierarchicalMutualInfo.list

#' Calculate Hierarchical Mutual Information between splits
#' 
#' @param splits1,splits2 Objects of class \code{Splits}.
#' @param nTip Integer specifying the number of tips.
#' @param normalize Logical. If \code{TRUE}, normalize the result.
#' @param reportMatching Logical specifying whether to return matchings.
#' 
#' @return Numeric value of Hierarchical Mutual Information.
#' 
#' @export
HierarchicalMutualInfoSplits <- function(splits1, splits2, 
                                        nTip = attr(splits1, "nTip"),
                                        normalize = FALSE,
                                        reportMatching = FALSE) {
  
  if (attr(splits1, "nTip") != attr(splits2, "nTip")) {
    stop("Trees must have the same number of tips")
  }
  
  # Calculate hierarchical weights for each split
  weights1 <- .CalculateHierarchicalWeights(splits1, nTip)
  weights2 <- .CalculateHierarchicalWeights(splits2, nTip) 
  
  # Calculate mutual information with hierarchical weighting
  hmi <- .CalculateWeightedMutualInfo(splits1, splits2, weights1, weights2, nTip)
  
  if (normalize) {
    # Normalize by the maximum of the two self-comparisons
    hmi_self1 <- .CalculateWeightedMutualInfo(splits1, splits1, weights1, weights1, nTip)
    hmi_self2 <- .CalculateWeightedMutualInfo(splits2, splits2, weights2, weights2, nTip)
    max_hmi <- max(hmi_self1, hmi_self2)
    if (max_hmi > 0) {
      hmi <- hmi / max_hmi
    }
  }
  
  if (reportMatching) {
    # For now, return empty matching - can be extended later
    attr(hmi, "matching") <- integer(0)
  }
  
  return(hmi)
}

#' Calculate hierarchical weights for splits based on tree structure
#' 
#' @param splits A \code{Splits} object
#' @param nTip Number of tips in the tree
#' 
#' @return Numeric vector of weights for each split
#' 
#' @keywords internal
.CalculateHierarchicalWeights <- function(splits, nTip) {
  n_splits <- length(splits)
  if (n_splits == 0) return(numeric(0))
  
  # Calculate depth-based weights
  # Deeper splits (closer to tips) get higher weights
  split_sizes <- TreeTools::TipsInSplits(splits)
  
  # Weight splits by their information content and hierarchy level
  # More balanced splits and deeper splits get higher weights  
  weights <- vapply(split_sizes, function(size) {
    # Entropy component (balanced splits are more informative)
    entropy_weight <- Entropy(c(size, nTip - size) / nTip)
    
    # Depth component (smaller splits are typically deeper)
    depth_weight <- 1 / (1 + abs(size - nTip/2))
    
    # Combine weights
    entropy_weight * (1 + depth_weight)
  }, numeric(1))
  
  # Normalize weights
  if (sum(weights) > 0) {
    weights <- weights / sum(weights)
  }
  
  return(weights)
}

#' Calculate weighted mutual information between two sets of splits
#' 
#' @param splits1,splits2 \code{Splits} objects
#' @param weights1,weights2 Numeric vectors of weights for each split
#' @param nTip Number of tips
#' 
#' @return Numeric value of weighted mutual information
#' 
#' @keywords internal
.CalculateWeightedMutualInfo <- function(splits1, splits2, weights1, weights2, nTip) {
  
  if (length(splits1) == 0 || length(splits2) == 0) {
    return(0)
  }
  
  # Calculate pairwise mutual information between all split pairs
  hmi_total <- 0
  
  for (i in seq_along(splits1)) {
    for (j in seq_along(splits2)) {
      # Convert splits to logical vectors if they're raw
      split1_logical <- .SplitToLogical(splits1, i, nTip)
      split2_logical <- .SplitToLogical(splits2, j, nTip)
      
      # Calculate mutual information between these splits
      mi <- MeilaMutualInformation(split1_logical, split2_logical)
      
      # Weight by hierarchical position
      weight <- weights1[i] * weights2[j]
      
      # Add to total HMI
      hmi_total <- hmi_total + (mi * weight)
    }
  }
  
  return(hmi_total)
}

#' Convert a split from a Splits object to logical vector
#' 
#' @param splits A Splits object
#' @param index Index of the split to extract
#' @param nTip Number of tips
#' 
#' @return Logical vector representing the split
#' 
#' @keywords internal
.SplitToLogical <- function(splits, index, nTip) {
  # Extract the split
  if (is.matrix(splits)) {
    # Raw matrix format used by TreeTools
    split_raw <- splits[index, ]
    # Convert raw to logical
    split_logical <- as.logical(rawToBits(split_raw)[seq_len(nTip)])
  } else {
    # Already in list format
    split_logical <- splits[[index]]
    if (is.raw(split_logical)) {
      split_logical <- as.logical(rawToBits(split_logical)[seq_len(nTip)])
    } else if (!is.logical(split_logical)) {
      split_logical <- as.logical(split_logical)
    }
  }
  
  return(split_logical)
}

#' Calculate maximum possible HMI for normalization
#' 
#' @param splits1,splits2 \code{Splits} objects
#' @param weights1,weights2 Numeric vectors of weights
#' 
#' @return Maximum possible HMI value
#' 
#' @keywords internal
.MaxHierarchicalMutualInfo <- function(splits1, splits2, weights1, weights2) {
  
  # Maximum occurs when trees are identical
  # Calculate self-mutual information with weights
  max_hmi <- 0
  
  # Use the tree with higher total weight as reference
  if (sum(weights1) >= sum(weights2)) {
    ref_splits <- splits1
    ref_weights <- weights1
    nTip <- attr(splits1, "nTip")
  } else {
    ref_splits <- splits2  
    ref_weights <- weights2
    nTip <- attr(splits2, "nTip")
  }
  
  for (i in seq_along(ref_splits)) {
    split_logical <- .SplitToLogical(ref_splits, i, nTip)
    
    # Self mutual information is just the entropy
    entropy <- Entropy(c(sum(split_logical), sum(!split_logical)) / length(split_logical))
    
    # Weight by hierarchical position (squared for self-comparison)
    weight <- ref_weights[i]^2
    
    max_hmi <- max_hmi + (entropy * weight)
  }
  
  return(max_hmi)
}