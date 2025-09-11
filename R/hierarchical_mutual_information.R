#' Hierarchical Mutual Information for phylogenetic trees
#'
#' Calculate the Hierarchical Mutual Information (\acronym{HMI})
#' between two phylogenetic trees, following the recursive algorithm from 
#' \insertCite{Perotti2015,Perotti2020;textual}{TreeDist}.
#' 
#' @details
#' This function implements the recursive Hierarchical Mutual Information algorithm
#' that considers the nested, hierarchical structure of phylogenetic trees when 
#' computing information measures. The algorithm converts trees to hierarchical 
#' partitions and computes mutual information recursively using natural logarithm.
#' 
#' The recursive \acronym{HMI} formula for internal nodes is:
#' \deqn{I(t,s) = ln(n_ts) - (H_us + H_tv - H_uv)/n_ts + mean(I_uv)}
#' 
#' Where:
#' \itemize{
#'   \item \eqn{n_ts} is the number of common elements between partitions
#'   \item \eqn{H_us, H_tv, H_uv} are entropy terms from child comparisons 
#'   \item \eqn{I_uv} is the recursive \acronym{HMI} for child pairs
#' }
#' 
#' @param tree1,tree2 Trees of class \code{phylo}, or lists of such trees.
#' If \code{tree2} is not provided, distances will be calculated between
#' each pair of trees in the list \code{tree1}.
#' @param normalize If `FALSE`, do not normalize the result.  If a function,
#' Normalize the result to range \[0,1\] by dividing by
#' `Func(SelfHMI(tree1), SelfHMI(tree2))`, where `Func()` = `max()` if 
#' `normalize == TRUE`, `normalize()` otherwise.
#' 
#' @return A numeric value representing the Hierarchical Mutual Information
#' between the input trees. Higher values indicate more shared 
#' hierarchical structure.
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
#' # HMI normalized against the mean information content of tree1 and tree2
#' HierarchicalMutualInfo(tree1, tree2, normalize = mean)
#' 
#' @references
#' \insertAllCited{}
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

#' @rdname HierarchicalMutualInfo
#' @export
HMI <- HierarchicalMutualInfo

#' @rdname HierarchicalMutualInfo
#' @export
HierarchicalMutualInformation <- HierarchicalMutualInfo

#' @return `SelfHMI()` returns the hierarchical mutual information of a tree
#' compared with itself, i.e. its hierarchical entropy (\acronym{HH}).
#' @examples 
#' # Normalized HMI above is equivalent to:
#' HMI(tree1, tree2) / mean(SelfHMI(tree1), SelfHMI(tree2))
#' @rdname HierarchicalMutualInfo
#' @export
SelfHMI <- function(tree) {
  part <- as.HPart(tree)
  HME_xptr(part)
}

#' @return `EHMI()` returns the expected \acronym{HMI} against a uniform
#' shuffling of element labels, estimated by performing Monte Carlo resampling
#' on the same hierarchical structure until the standard error of the
#' estimate falls below `tolerance`.
#' The attributes of the returned object list the variance (`var`),
#' standard deviation (`sd`), standard error of the mean (`sem`) and 
#' relative error (`relativeError`) of the estimate, and the number of Monte
#' Carlo samples used to obtain it (`samples`).
#' @examples
#' # Expected mutual info for this pair of hierarchies
#' EHMI(tree1, tree2, tolerance = 0.1)
#' @rdname HierarchicalMutualInfo
#' @export
EHMI <- function(tree1, tree2, tolerance = 0.01, minResample = 36) {
  EHMI_xptr(as.HPart(tree1), as.HPart(tree2), as.numeric(tolerance),
                as.integer(minResample))
}

.AHMISEM <- function(hmi, M, ehmi, ehmi_sem) {
  deriv <- (hmi - M) / (M - ehmi)^2
  abs(deriv) * ehmi_sem
}

#' @return `AHMI()` returns the adjusted \acronym{HMI}, normalized such that
#' zero corresponds to the expected \acronym{HMI} given a random shuffling
#' of elements on the same hierarchical structure.  The attribute `sem` gives
#' the standard error of the estimate.
#' @examples
#' # The adjusted HMI normalizes against this expectation
#' AHMI(tree1, tree2, tolerance = 0.1)
#' @rdname HierarchicalMutualInfo
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
