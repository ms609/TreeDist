#' Hierarchical Mutual Information for phylogenetic trees
#'
#' Calculate the Hierarchical Mutual Information (\acronym{HMI})
#' between two phylogenetic trees, following the recursive algorithm of
#' \insertCite{Perotti2020;textual}{TreeDist}.
#' 
#' @details
#' `HierarchicalMutualInfo()` computes the hierarchical mutual content of trees
#' \insertCite{Perotti2015,Perotti2020}{TreeDist}, which accounts for the
#' non-independence of information represented by nested splits.
#' 
#' `tree` is converted to a set of hierarchical partitions, and the mutual
#' information (in bits) is computed recursively; the contribution of a node is
#' given by:
#' 
#' \deqn{I(t,s) = \log_2(n_{ts}) - \dfrac{H_{us} + H_{tv} - H_{uv}}{n_{ts}} +
#' \text{mean}(I_{uv})}
#' 
#' Where:
#' \itemize{
#'   \item \eqn{n_{ts}} is the number of common elements between partitions
#'   \item \eqn{H_{us}, H_{tv}, H_{uv}} are entropy terms from child comparisons 
#'   \item \eqn{I_{uv}} is the recursive \acronym{HMI} for child pairs
#' }
#' 
#' @param tree,tree1,tree2 An object that can be coerced to an [`HPart`] 
#' object, or (soon) a list of such objects.
#' (Not yet implemented: ) If \code{tree2} is not provided, distances will be
#' calculated between each pair of trees in the list \code{tree1}.
#' @param normalize If `FALSE`, return the raw \acronym{HMI}, in bits.
#' If `TRUE`, normalize to range \[0,1\] by dividing by
#' `max(SelfHMI(tree1), SelfHMI(tree2))`.
#' If a function, divide by `normalize(SelfHMI(tree1), SelfHMI(tree2))`.
#' 
#' @return `HierarchicalMutualInfo()` returns a numeric value representing the
#' Hierarchical Mutual Information between the input trees, in bits,
#' normalized as specified.
#' Higher values indicate more shared hierarchical structure.
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
      hmi / log(2)
    } else {
      if (isTRUE(normalize)) {
        normalize <- max
      }
      if (!is.function(normalize)) {
        stop("`normalize` must be logical, or a function")
      }
      denom <- normalize(HH_xptr(hp1), HH_xptr(hp2))
      hmi / denom
    }
  }
}

#' @keywords internal
#' @export
HierarchicalMutualInformation <- HierarchicalMutualInfo

#' @rdname HierarchicalMutualInfo
#' @export
HMI <- HierarchicalMutualInfo

#' @return `SelfHMI()` returns the hierarchical mutual information of a tree
#' compared with itself, i.e. its hierarchical entropy (\acronym{HH}).
#' @examples 
#' # Normalized HMI above is equivalent to:
#' HMI(tree1, tree2) / mean(SelfHMI(tree1), SelfHMI(tree2))
#' @rdname HierarchicalMutualInfo
#' @export
SelfHMI <- function(tree) {
  part <- as.HPart(tree)
  HH_xptr(part) / log(2)
}

CharH <- function(tree) {
  part <- as.HPart(tree)
  H_xptr(part) / log(2)
}

#' @export
#' @keywords internal
HH <- SelfHMI

#' @param precision Numeric; Monte Carlo sampling will terminate once the
#' relative standard error falls below this value.
#' @param minResample Integer specifying minimum number of Monte Carlo samples
#' to conduct.  Avoids early termination when sample size is too small to
#' reliably estimate the standard error of the mean.
#' @return `EHMI()` returns the expected \acronym{HMI} against a uniform
#' shuffling of element labels, estimated by performing Monte Carlo resampling
#' on the same hierarchical structure until the relative standard error of the
#' estimate falls below `precision`.
#' The attributes of the returned object list the variance (`var`),
#' standard deviation (`sd`), standard error of the mean (`sem`) and 
#' relative error (`relativeError`) of the estimate, and the number of Monte
#' Carlo samples used to obtain it (`samples`).
#' @examples
#' # Expected mutual info for this pair of hierarchies
#' EHMI(tree1, tree2, precision = 0.1)
#' @rdname HierarchicalMutualInfo
#' @export
EHMI <- function(tree1, tree2, precision = 0.01, minResample = 36) {
  EHMI_xptr(as.HPart(tree1), as.HPart(tree2), as.numeric(precision),
                as.integer(minResample)) / log(2)
}

.AHMISEM <- function(hmi, M, ehmi, ehmi_sem) {
  deriv <- (hmi - M) / (M - ehmi)^2
  abs(deriv) * ehmi_sem
}

#' @details `AHMI()` calculates the adjusted hierarchical mutual information:
#' \deqn{\text{AHMI}(t, s) = \dfrac{I(t, s) - \hat{I}(t, s)}{
#'  \text{mean}(H(t), H(s)) - \hat{I}(t, s)}}
#' Where:
#' - \eqn{I(t, s)} is the hierarchical mutual information between `tree1` and
#'   `tree2`
#' - \eqn{\hat{I}(t, s)} is the expected \acronym{HMI} between `tree1` and
#'   `tree2`, estimated by Monte Carlo sampling
#' - \eqn{H(t), H(s)} is the entropy (self-mutual information) of each tree
#' @param Mean Function by which to combine the self-information of the 
#' two input hierarchies, in order to normalize the \acronym{HMI}.
#' @return `AHMI()` returns the adjusted \acronym{HMI}, normalized such that
#' zero corresponds to the expected \acronym{HMI} given a random shuffling
#' of elements on the same hierarchical structure.  The attribute `sem` gives
#' the standard error of the estimate.
#' @examples
#' # The adjusted HMI normalizes against this expectation
#' AHMI(tree1, tree2, precision = 0.1)
#' @rdname HierarchicalMutualInfo
#' @export
AHMI <- function(tree1, tree2, Mean = max, precision = 0.01, minResample = 36) {
  hp1 <- as.HPart(tree1)
  hp2 <- as.HPart(tree2, hp1)
  
  ehmi <- EHMI_xptr(hp1, hp2, as.numeric(precision), as.integer(minResample))
  hmi <- HMI_xptr(hp1, hp2)
  hh1 <- HH_xptr(hp1)
  hh2 <- HH_xptr(hp2)
  M <- Mean(hh1, hh2)
  
  num <- hmi - ehmi[[1]]
  denom <- M - ehmi[[1]]
  # Return:
  structure(if (num < sqrt(.Machine$double.eps)) 0 else num / denom,
            sem = .AHMISEM(hmi, M, ehmi[[1]], attr(ehmi, "sem")))
}
