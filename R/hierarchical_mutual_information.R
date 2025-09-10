#' Hierarchical Mutual Information for phylogenetic trees
#'
#' Calculate the Hierarchical Mutual Information (HMI) between two phylogenetic
#' trees, following the recursive algorithm from Perotti et al. (2015).
#' 
#' @details
#' This function implements the recursive Hierarchical Mutual Information algorithm
#' that considers the nested, hierarchical structure of phylogenetic trees when 
#' computing information measures. The algorithm converts trees to hierarchical 
#' partitions and computes mutual information recursively using natural logarithm.
#' 
#' The recursive HMI formula for internal nodes is:
#' I(t,s) = ln(n_ts) - (H_us + H_tv - H_uv)/n_ts + mean(I_uv)
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
#' @param normalize If `FALSE`, do not normalize the result.  If a function,
#' Normalize the result to range [0,1] by dividing by
#' `Func(SelfHMI(tree1), SelfHMI(tree2))`, where `Func()` = `max()` if 
#' `normalize == TRUE`, `normalize()` otherwise.
#' @param reportMatching Logical specifying whether to return the clade
#' matchings as an attribute of the score.
#' 
#' @return A numeric value representing the Hierarchical Mutual Information
#' between the input trees. Higher values indicate more shared 
#' hierarchical structure.
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
#' HierarchicalMutualInfo(bal6, pec6)  # Returns approximately 0.22
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
HierarchicalMutualInfo <- function(tree1, tree2 = NULL, normalize = FALSE) {
  hp1 <- as.HPart(tree1)
  if (is.null(tree2)) {
    if (isFALSE(normalize)) {
      SelfHMI(hp1)
    } else {
      warning("Normalized self-information == 1; did you mean to provide tree2?")
      1
    }
  } else {
    hp2 <- as.HPart(tree2, tree1)
    hmi <- HMI_xptr(hp1, hp2)
    if (isFALSE(normalize)) {
      hmi
    } else {
      if (isTRUE(normalize)) {
        normalize <- max
      }
      if (!is.function(normalize)) {
        stop("`normalize` must be logical, or a function")
      }
      denom <- normalize(SelfHMI(hp1), SelfHMI(hp2))
      hmi / denom
    }
  }
}

#' @export
HMI <- HierarchicalMutualInfo

#' @export
SelfHMI <- function(tree) {
  part <- as.HPart(tree)
  HMI(part, part)[[2]]
}

#' @export
SelfHMI <- function(tree) {
  part <- as.HPart(tree)
  HME_xptr(part)
}

#' @export
EHMI <- function(tree1, tree2, tolerance = 0.01, minResample = 36) {
  EHMI_xptr(as.HPart(tree1), as.HPart(tree2), as.numeric(tolerance),
                as.integer(minResample))
}

.AHMISEM <- function(hmi, M, ehmi, ehmi_sem) {
  deriv <- (hmi - M) / (M - ehmi)^2
  abs(deriv) * ehmi_sem
}

#' @export
AHMI <- function(tree1, tree2, Mean = max, tolerance = 0.01, minResample = 36) {
  hp1 <- as.HPart(tree1)
  hp2 <- as.HPart(tree2, hp1)
  
  ehmi <- EHMI_xptr(hp1, hp2, as.numeric(tolerance), as.integer(minResample))
  hmi <- HMI_xptr(hp1, hp2)
  hh1 <- HME_xptr(hp1)
  hh2 <- HME_xptr(hp2)
  M <- Mean(hh1, hh2)
  
  # Return:
  structure((hmi - ehmi[[1]]) / (M - ehmi[[1]]),
            sem = .AHMISEM(hmi, M, ehmi[[1]], attr(ehmi, "sem")))
}
