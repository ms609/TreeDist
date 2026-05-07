#' @rdname TreeDistance
#' @export
MatchingSplitInfo <- function(tree1, tree2 = NULL, normalize = FALSE,
                               reportMatching = FALSE, diag = TRUE) {
  unnormalized <- CalculateTreeDistance(MatchingSplitInfoSplits, tree1,
                                        tree2, reportMatching)
  
  if (diag && is.null(tree2)) {
    unnormalized <- as.matrix(unnormalized)
    diag(unnormalized) <- SplitwiseInfo(tree1)
    tree2 <- tree1
  }
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = SplitwiseInfo, Combine = .PairMean)
}

#' @rdname TreeDistance
#' @export
MatchingSplitInfoDistance <- function(tree1, tree2 = NULL, 
                                       normalize = FALSE,
                                       reportMatching = FALSE) {
  
  # Fast path (all-pairs): same tips, no matching — avoids duplicate as.Splits()
  fast <- .FastDistPath(tree1, tree2, reportMatching,
                        cpp_msi_all_pairs,
                        cpp_splitwise_info_batch)
  if (!is.null(fast)) {
    msi <- fast[["info"]]
    treesIndependentInfo <- .PairwiseSums(fast[["entropies"]])
    
    ret <- .FloorNumericalNoise(treesIndependentInfo - msi - msi, treesIndependentInfo)
    ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                         infoInBoth = treesIndependentInfo,
                         InfoInTree = SplitwiseInfo, Combine = "+")
    attributes(ret) <- attributes(msi)
    return(ret)
  }

  # Fast path (cross-pairs): same tips, no matching — avoids duplicate as.Splits()
  fast_many <- .FastManyManyPath(tree1, tree2, reportMatching,
                                 cpp_msi_cross_pairs,
                                 cpp_splitwise_info_batch)
  if (!is.null(fast_many)) {
    msi <- fast_many[["dists"]]
    info1 <- fast_many[["info1"]]
    info2 <- fast_many[["info2"]]
    treesIndependentInfo <- outer(info1, info2, "+")

    ret <- .FloorNumericalNoise(treesIndependentInfo - msi - msi, treesIndependentInfo)
    ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                         infoInBoth = treesIndependentInfo,
                         InfoInTree = SplitwiseInfo, Combine = "+")
    return(ret)
  }

  msi <- MatchingSplitInfo(tree1, tree2, normalize = FALSE, diag = FALSE,
                           reportMatching = reportMatching)

  treesIndependentInfo <- .MaxValue(tree1, tree2, SplitwiseInfo)
  ret <- .FloorNumericalNoise(treesIndependentInfo - msi - msi, treesIndependentInfo)
  ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                       infoInBoth = treesIndependentInfo,
                       InfoInTree = SplitwiseInfo, Combine = "+")
  attributes(ret) <- attributes(msi)
  # Return:
  ret
}


#' @rdname TreeDistance
#' @export
MatchingSplitInfoSplits <- function(splits1, splits2,
                                     nTip = attr(splits1, "nTip"),
                                     reportMatching = FALSE) {
  
  GeneralizedRF(splits1, splits2, nTip, cpp_msi_distance, maximize = TRUE,
                reportMatching = reportMatching)
}
