#' Robinson&ndash;Foulds distances, with adjustments for phylogenetic information
#' content
#' 
#' `RobinsonFoulds()` calculates the Robinson&ndash;Foulds distance
#' \insertCite{Robinson1981}{TreeDist}, or the corresponding similarity measure.
#' `InfoRobinsonFoulds()` weights splits according to their phylogenetic
#' information content \insertCite{@§2.1 in @SmithDist}{TreeDist}.
#' Optionally, the matching between identical splits may reported.
#' Generalized Robinson&ndash;Foulds distances (see [`TreeDistance()`])
#' are better suited to most use cases
#' \insertCite{SmithDist,SmithSpace}{TreeDist}.
#' 
#' `RobinsonFoulds()` calculates the standard Robinson&ndash;Foulds distance,
#' i.e. the number of splits that occur in one tree but not the other.
#' `InfoRobinsonFoulds()` calculates the tree similarity or distance by summing 
#' the phylogenetic information content of all splits that are (or are not)
#' identical in both trees.  Consequently, splits that are more likely
#' to be identical by chance alone make a smaller contribution to overall
#' tree distance, because their similarity is less remarkable.
#' 
#' Rapid comparison between multiple pairs of trees employs the
#' \insertCite{Day1985;textual}{TreeDist} linear-time algorithm.
#' 
#' @inheritParams TreeDistance
#' @param similarity Logical specifying whether to report the result as a tree
#' similarity, rather than a difference.
#' 
#' @templateVar returns `RobinsonFoulds()` and `InfoRobinsonFoulds()` return
#' @template distReturn
#' @return If `reportMatching = TRUE`, the `pairScores` attribute 
#' returns a logical matrix specifying whether each pair of splits is identical.
#' 
#' 
#' @section Normalization:
#' 
#' - `RobinsonFoulds()` is normalized against the total number of splits that
#'  are present.
#'  
#' - `InfoRobinsonFoulds()` is normalized against the sum of the phylogenetic 
#' information of all splits in each tree, treated independently.
#'  
#' @references \insertAllCited{}
#' 
#' @examples
#'  # For BalancedTree, PectinateTree, as.phylo:
#' library("TreeTools", quietly = TRUE)
#' balanced7 <- BalancedTree(7)
#' pectinate7 <- PectinateTree(7)
#' RobinsonFoulds(balanced7, pectinate7)
#' RobinsonFoulds(balanced7, pectinate7, normalize = TRUE)
#' VisualizeMatching(RobinsonFouldsMatching, balanced7, pectinate7)
#' 
#' InfoRobinsonFoulds(balanced7, pectinate7)
#' VisualizeMatching(InfoRobinsonFoulds, balanced7, pectinate7)
#' @template MRS
#' 
#' 
#' @family tree distances
#' @seealso Display paired splits: [`VisualizeMatching()`]
#' 
#' @export
#' @encoding UTF-8
#' @name Robinson-Foulds

#' @aliases RobinsonFouldsInfo
#' @rdname Robinson-Foulds
InfoRobinsonFoulds <- function(tree1, tree2 = NULL, similarity = FALSE,
                                normalize = FALSE, reportMatching = FALSE) {
  
  # Fast path for distance (not similarity): avoids duplicate as.Splits()
  if (!similarity) {
    # All-pairs fast path
    fast <- .FastDistPath(tree1, tree2, reportMatching,
                          cpp_rf_info_all_pairs,
                          cpp_splitwise_info_batch)
    if (!is.null(fast)) {
      treesIndependentInfo <- .PairwiseSums(fast[["entropies"]])
      unnormalized <- treesIndependentInfo - fast[["info"]] - fast[["info"]]
      unnormalized[unnormalized < .Machine[["double.eps"]] ^ 0.5] <- 0
      ret <- NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                           InfoInTree = SplitwiseInfo, Combine = "+")
      attributes(ret) <- attributes(fast[["info"]])
      return(ret)
    }
    
    # Cross-pairs fast path
    fast_many <- .FastManyManyPath(tree1, tree2, reportMatching,
                                   cpp_rf_info_cross_pairs,
                                   cpp_splitwise_info_batch)
    if (!is.null(fast_many)) {
      irf <- fast_many[["dists"]]
      info1 <- fast_many[["info1"]]
      info2 <- fast_many[["info2"]]
      treesIndependentInfo <- outer(info1, info2, "+")
      
      unnormalized <- treesIndependentInfo - irf - irf
      unnormalized[unnormalized < .Machine[["double.eps"]] ^ 0.5] <- 0
      ret <- NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                           InfoInTree = SplitwiseInfo, Combine = "+")
      return(ret)
    }
  }
  
  unnormalized <- CalculateTreeDistance(InfoRobinsonFouldsSplits, tree1, tree2, 
                                        reportMatching) * 2
  
  if (!similarity) {
    unnormalized <- .MaxValue(tree1, tree2, SplitwiseInfo) - unnormalized
  }
  
  # In case of floating point inaccuracy
  unnormalized[unnormalized < .Machine[["double.eps"]] ^ 0.5] <- 0
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = SplitwiseInfo, Combine = "+")
}

#' @export
RobinsonFouldsInfo <- InfoRobinsonFoulds

#' @rdname Robinson-Foulds
#' @inheritParams SharedPhylogeneticInfoSplits
#' @export
InfoRobinsonFouldsSplits <- function(splits1, splits2, 
                                      nTip = attr(splits1, "nTip"),
                                      reportMatching = FALSE) {
  
  GeneralizedRF(splits1, splits2, nTip, cpp_robinson_foulds_info,
                maximize = FALSE, reportMatching = reportMatching)
}

#' @rdname Robinson-Foulds
#' @importFrom TreeTools as.ClusterTable NSplits
#' @export
RobinsonFoulds <- function(tree1, tree2 = NULL, similarity = FALSE,
                            normalize = FALSE, reportMatching = FALSE) {
  if (is.null(tree2)) {
    ct <- as.ClusterTable(tree1)
    rf <- robinson_foulds_all_pairs(if (is.list(ct)) ct else list(ct))
    if (similarity) {
      unnormalized <- structure(rf + rf, Size = length(tree1), class = "dist")
    } else {
      splits <- NSplits(tree1)
      nSplits <- outer(splits, splits, "+")
      unnormalized <- structure(nSplits[lower.tri(nSplits)] - rf - rf,
                                Size = length(tree1),
                                class = "dist")
    }
  } else {
    # Fast cross-pairs path: batch C++ via ClusterTable (Day 1985).
    # Only applicable when both inputs are lists of trees with matching
    # tip labels and reportMatching is not requested.
    fast <- NULL
    if (!reportMatching &&
        !inherits(tree1, c("phylo", "Splits")) &&
        !inherits(tree2, c("phylo", "Splits")) &&
        is.null(getOption("TreeDist-cluster"))) {
      lab1 <- TipLabels(tree1)
      lab2 <- TipLabels(tree2)
      if (is.list(lab1)) lab1 <- lab1[[1]]
      if (is.list(lab2)) lab2 <- lab2[[1]]
      if (setequal(lab1, lab2)) {
        ct1 <- as.ClusterTable(tree1, tipLabels = lab1)
        ct2 <- as.ClusterTable(tree2, tipLabels = lab1)
        shared <- robinson_foulds_cross_pairs(
          if (is.list(ct1)) ct1 else list(ct1),
          if (is.list(ct2)) ct2 else list(ct2)
        )
        splits1 <- NSplits(tree1)
        splits2 <- NSplits(tree2)
        if (similarity) {
          fast <- shared + shared
        } else {
          fast <- outer(splits1, splits2, "+") - shared - shared
        }
        dimnames(fast) <- list(names(tree1), names(tree2))
      }
    }
    if (is.null(fast)) {
      unnormalized <- CalculateTreeDistance(RobinsonFouldsSplits, tree1, tree2,
                                            reportMatching)
      if (similarity) {
        unnormalized <- .MaxValue(tree1, tree2, NSplits) - unnormalized
      }
    } else {
      unnormalized <- fast
    }
  }
  
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = NSplits, Combine = `+`)
}

#' @describeIn Robinson-Foulds Matched splits, intended for use with 
#' [`VisualizeMatching()`].
#' @param \dots Not used.
#' @importFrom TreeTools NSplits
#' @export
RobinsonFouldsMatching <- function(tree1, tree2, similarity = FALSE,
                                    normalize = FALSE, ...) {
  ret <- CalculateTreeDistance(RobinsonFouldsSplits, tree1, tree2,
                               reportMatching = TRUE)

  ret <- .MaxValue(tree1, tree2, NSplits) - ret
  
  attr(ret, "pairScores") <- !attr(ret, "pairScores")
  
  # Return:
  ret
}

#' @rdname Robinson-Foulds
#' @inheritParams SharedPhylogeneticInfoSplits
#' @export
RobinsonFouldsSplits <- function(splits1, splits2,
                                  nTip = attr(splits1, "nTip"),
                                  reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, nTip, cpp_robinson_foulds_distance,
                maximize = FALSE, reportMatching = reportMatching)
}
