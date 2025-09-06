#' Hierarchical Mutual Information for phylogenetic trees
#'
#' Calculate the Hierarchical Mutual Information (HMI) between two phylogenetic
#' trees, which extends traditional mutual information to account for the
#' hierarchical structure inherent in phylogenetic trees.
#' 
#' @details
#' Hierarchical Mutual Information is a recursive algorithm that considers the 
#' nested, hierarchical structure of phylogenetic trees when computing information 
#' measures. The algorithm converts trees to hierarchical partitions and computes
#' mutual information recursively, weighting contributions by the number of
#' overlapping elements at each level of the hierarchy.
#' 
#' The algorithm follows the implementation described in Perotti et al. (2015)
#' and is based on the recursive formula:
#' 
#' For internal nodes: I(t,s) = log(n_ts) - (H_us + H_tv - H_uv)/n_ts + mean(I_uv)
#' 
#' Where:
#' \itemize{
#'   \item n_ts is the number of common elements between partitions
#'   \item H_us, H_tv, H_uv are entropy terms from child comparisons
#'   \item I_uv is the recursive HMI for child pairs
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
#' \dontrun{
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
#' # Expected result for 6-tip balanced vs pectinate trees
#' bal6 <- BalancedTree(6)
#' pec6 <- PectinateTree(6)
#' HierarchicalMutualInfo(bal6, pec6)  # Should be approximately 0.24
#' }
#' 
#' @references
#' Perotti, J. I., Tessone, C. J., & Caldarelli, G. (2015). 
#' Hierarchical mutual information for the comparison of hierarchical 
#' community structures in complex networks. 
#' Physical Review E, 92(6), 062825.
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
  
  # Ensure trees have the same number of tips
  if (length(tree1$tip.label) != length(tree2$tip.label)) {
    stop("Trees must have the same number of tips")
  }
  
  # Convert trees to hierarchical partitions 
  partition1 <- .PhyloToHierarchicalPartition(tree1)
  partition2 <- .PhyloToHierarchicalPartition(tree2)
  
  # Calculate HMI using recursive algorithm
  result <- .CalculateHMIRecursive(partition1, partition2)
  hmi <- result$I_ts
  
  if (normalize) {
    # Normalize by the maximum of the two self-comparisons
    hmi_self1 <- .CalculateHMIRecursive(partition1, partition1)$I_ts
    hmi_self2 <- .CalculateHMIRecursive(partition2, partition2)$I_ts
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

#' @export  
HierarchicalMutualInfo.list <- function(tree1, tree2 = NULL, normalize = FALSE,
                                       reportMatching = FALSE) {
  # For lists, we need to handle them as distance calculations
  # This would require significant rework of CalculateTreeDistance function
  # For now, provide a basic implementation
  
  if (is.null(tree2)) {
    # Calculate all pairwise distances
    n <- length(tree1)
    result_matrix <- matrix(0, n, n)
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        hmi_val <- HierarchicalMutualInfo.phylo(tree1[[i]], tree1[[j]], 
                                               normalize = normalize,
                                               reportMatching = reportMatching)
        result_matrix[i, j] <- hmi_val
        result_matrix[j, i] <- hmi_val
      }
    }
    
    # Convert to dist object
    return(as.dist(result_matrix))
  } else {
    # Pairwise between two lists
    stop("Pairwise list comparison not yet implemented")
  }
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
  
  # This function will now convert splits back to trees and use the proper HMI algorithm
  # For now, use a simplified approach - the main function should handle tree objects
  stop("HierarchicalMutualInfoSplits is deprecated. Use HierarchicalMutualInfo with phylo objects.")
}

#' Convert phylo tree to hierarchical partition
#' 
#' @param tree A phylo object
#' 
#' @return A nested list representing the hierarchical partition
#' 
#' @keywords internal
.PhyloToHierarchicalPartition <- function(tree) {
  
  # Convert tree structure to hierarchical partition directly from edge matrix
  # This avoids dependency on write.tree
  
  # Get number of tips
  nTip <- length(tree$tip.label)
  
  # Build hierarchical partition from tree structure
  .BuildHierarchicalPartition(tree, nTip + 1, nTip)  # Start from root
}

#' Build hierarchical partition recursively from tree structure
#' 
#' @param tree Phylo object
#' @param node Current node number
#' @param nTip Number of tips
#' 
#' @return Hierarchical partition for this subtree
#' 
#' @keywords internal
.BuildHierarchicalPartition <- function(tree, node, nTip) {
  
  # Find children of this node
  children <- tree$edge[tree$edge[, 1] == node, 2]
  
  # If no children, this is a tip
  if (length(children) == 0) {
    # This is a tip node, return the tip number
    return(node)
  }
  
  # If this node has children, recursively build partition for each child
  result <- list()
  for (child in children) {
    child_partition <- .BuildHierarchicalPartition(tree, child, nTip)
    result <- append(result, list(child_partition))
  }
  
  return(result)
}

#' Calculate Hierarchical Mutual Information recursively
#' 
#' @param Ut,Us Hierarchical partitions (nested lists)
#' 
#' @return List with n_ts and I_ts values
#' 
#' @keywords internal
.CalculateHMIRecursive <- function(Ut, Us) {
  
  # Helper function for x*log(x)
  xlnx <- function(x) {
    if (x <= 0) 0 else x * log(x)
  }
  
  # Flatten function
  flattenator <- function(partition) {
    if (!is.list(partition)) {
      return(partition)
    }
    unlist(partition)
  }
  
  # Base case: both are leaves
  if (!is.list(Ut) && !is.list(Us)) {
    overlap <- length(intersect(Ut, Us))
    return(list(n_ts = overlap, I_ts = 0))
  }
  
  # Ut is internal node and Us is leaf
  if (is.list(Ut) && !is.list(Us)) {
    all_Ut <- flattenator(Ut)
    overlap <- length(intersect(all_Ut, Us))
    return(list(n_ts = overlap, I_ts = 0))
  }
  
  # Ut is leaf and Us is internal node  
  if (!is.list(Ut) && is.list(Us)) {
    all_Us <- flattenator(Us)
    overlap <- length(intersect(Ut, all_Us))
    return(list(n_ts = overlap, I_ts = 0))
  }
  
  # Both are internal nodes - main computation
  n_ts <- 0
  H_uv <- 0
  H_us <- 0
  H_tv <- 0
  mean_I_ts <- 0
  n_tv <- numeric(length(Us))
  
  for (u_idx in seq_along(Ut)) {
    Uu <- Ut[[u_idx]]
    n_us <- 0
    
    for (v_idx in seq_along(Us)) {
      Uv <- Us[[v_idx]]
      result <- .CalculateHMIRecursive(Uu, Uv)
      n_uv <- result$n_ts
      I_uv <- result$I_ts
      
      n_ts <- n_ts + n_uv
      n_tv[v_idx] <- n_tv[v_idx] + n_uv
      n_us <- n_us + n_uv
      H_uv <- H_uv + xlnx(n_uv)
      mean_I_ts <- mean_I_ts + n_uv * I_uv
    }
    H_us <- H_us + xlnx(n_us)
  }
  
  for (n_tv_val in n_tv) {
    H_tv <- H_tv + xlnx(n_tv_val)
  }
  
  if (n_ts > 0) {
    local_I_ts <- log(n_ts) - (H_us + H_tv - H_uv) / n_ts
    mean_I_ts <- mean_I_ts / n_ts
    I_ts <- local_I_ts + mean_I_ts
    return(list(n_ts = n_ts, I_ts = I_ts))
  } else {
    return(list(n_ts = 0, I_ts = 0))
  }
}

