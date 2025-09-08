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
#' @param normalize Logical. If \code{TRUE}, normalize the result to range [0,1].
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
  UseMethod("HierarchicalMutualInfo")
}

XLnX <- function(x) {
  ifelse(x > 0, x * log(x), 0)
}

#' @export
ReplicateHPart <- function(x, d) {
  rapply(x, function(x) d[[x]], how = "replace")
}

#' @importFrom stats setNames
#' @export
ShuffleHPart <- function(x) {
  labels <- as.character(TipLabels(x))
  d <- setNames(sample(labels), labels)
  ReplicateHPart(x, d)
}


#'   Computes the hierarchical mutual information between two hierarchical partitions.
#' @return Returns 
#' n_ts,HMI(Ut,Us) : where n_ts is the number of common elements between the hierarchical #' partitions Ut and Us.
#'
#' NOTE: We label by u,v the children of t,s respectively.
#' @export
HMI <- function(Ut, Us) {
  if (is.list(Ut[[1]])) {
    if (is.list(Us[[1]])) {
      # Ut and Us are both internal nodes since they contain other lists.
      n_ts = 0
      H_uv = 0
      H_us = 0
      H_tv = 0
      n_tv = integer(length(Us))
      mean_I_ts = 0
      for (Uu in Ut) {
        n_us = 0
        for (v in seq_along(Us)) {
          Uv <- Us[[v]]
          niUV <- HMI(Uu, Uv)
          n_uv <- niUV[[1]]
          I_uv <- niUV[[2]]
          n_ts <- n_ts + n_uv
          n_tv[[v]] <- n_tv[[v]] + n_uv
          n_us <- n_us + n_uv
          H_uv <- H_uv + XLnX(n_uv)
          mean_I_ts <- mean_I_ts + (n_uv * I_uv)
        }
        H_us <- H_us + XLnX(n_us)
      }
      for (.n_tv in n_tv) {
        H_tv <- H_tv + XLnX(.n_tv)
      }
      if (n_ts > 0) {
        local_I_ts <- log(n_ts) - (H_us + H_tv - H_uv) / n_ts
        mean_I_ts <- mean_I_ts / n_ts
        I_ts <- local_I_ts + mean_I_ts
        c(n_ts, I_ts)
      } else {
        c(0, 0)
      }
    } else {
      # Ut is internal node and Us is leaf
      c(length(intersect(unlist(Ut, recursive = TRUE), Us)), 0)
    }
  } else {
      if (is.list(Us)) {
        # Ut is leaf and Us internal node
        c(length(intersect(unlist(Us, recursive = TRUE), Ut)), 0)
      } else {
        # Both Ut and Us are leaves
        c(length(intersect(Ut, Us)), 0)
      }
  }
}

# TODO implement more efficiently
#' @export
SelfHMI <- function(tree) {
  part <- as.HPart(tree)
  HMI(part, part)[[2]]
}

#' @export
NHMI <- function(tree1, tree2) {
  part1 <- as.HPart(tree1)
  part2 <- as.HPart(tree2)
  gm <- mean(SelfHMI(part1), SelfHMI(part2))
  if (gm > 0) {
    HMI(part1, part2)[[2]] / gm
  } else {
    0
  }
}

EHMI <- function(tree1, tree2, tolerance = 0.01, minResample = 36) {
  if (minResample < 2) {
    stop("Must perform at least one resampling")
  }
  
  part1 <- as.HPart(tree1)
  part2 <- as.HPart(tree2)
  
  part1 <- rapply(part1, as.character, how = "replace")
  part2 <- rapply(part2, as.character, how = "replace")
  
  relativeError <- 2 * tolerance
  
  runMean <- 0
  runS <- 0
  runN <- 0
  
  progBar <- cli::cli_progress_bar("Sampling", total = NA, format = "{cli::pb_spin} Sample {runN}: {signif(runMean, 3)} ± {signif(runSEM, 3)} ({signif(relativeError * 100, 3)}%)")

  while(relativeError > tolerance || runN < minResample) {
    shuf1 <- ShuffleHPart(part1)
    x <- HMI(shuf1, part2)[[2]]
    
    runN <- runN + 1
    oldMean <- runMean
    runMean <- runMean + (x - runMean) / runN
    runS <- runS + (x - oldMean) * (x - runMean)
    runVar <- runS / (runN - 1)
    runSD <- sqrt(runVar)
    runSEM <- runSD / sqrt(runN)
    tolSD <- 0.05
    relativeError <- runSEM / (abs(runMean) + tolSD)
    cli::cli_progress_update(id = progBar,
                             status = list(runN = runN, runMean = runMean,
                                           runSEM = runSEM,
                                           relativeError = relativeError))
  }
  cli::cli_progress_done()
  
  structure(runMean, var = runVar, sd = runSD, sem = runSEM,
            relativeError = relativeError)
}

#' @export
AHMI <- function(tree1, tree2, Mean = max, tolerance = 0.01, minResample = 36) {
  hp1 <- as.HPart(tree1)
  hp2 <- as.HPart(tree2)
  ehmi <- EHMI(hp1, hp2, tolerance = tolerance, minResample = minResample)[[1]]
  # Return:
  (HMI(hp1, hp2)[[2]] - ehmi) / (Mean(SelfHMI(hp1), SelfHMI(hp2)) - ehmi)
}
