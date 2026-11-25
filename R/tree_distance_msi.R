#' @rdname TreeDistance
#' @export
MatchingSplitInfo <- function(tree1, tree2 = NULL, normalize = FALSE,
                               reportMatching = FALSE, diag = TRUE) {
  fast <- .FastSimilarity(tree1, tree2, normalize, reportMatching, diag,
                          cpp_msi_all_pairs, cpp_splitwise_info_batch,
                          cpp_msi_cross_pairs)
  if (!is.null(fast)) {
    return(NormalizeInfo(fast[["scores"]], tree1, tree2, SplitwiseInfo,
                         infoInBoth = fast[["infoInBoth"]], how = normalize,
                         Combine = .PairMean))
  }
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
    info <- fast[["entropies"]]
    treesIndependentInfo <- .PairwiseSums(info)
    normInfo <- .DistInfoToNormalize(info, info, normalize,
                                     treesIndependentInfo, lower = TRUE)
    
    ret <- .FloorNumericalNoise(treesIndependentInfo - msi - msi, treesIndependentInfo)
    ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                         infoInBoth = normInfo,
                         InfoInTree = SplitwiseInfo, Combine = "+")
    attributes(ret) <- attributes(msi)
    return(ret)
  }

  # Fast path (cross-pairs): same tips, no matching — avoids duplicate as.Splits()
  fast_many <- .FastManyManyPath(tree1, tree2, reportMatching,
                                 cpp_msi_cross_pairs)
  if (!is.null(fast_many)) {
    msi <- fast_many[["dists"]]
    treesIndependentInfo <- fast_many[["infoBoth"]]
    normInfo <- .DistInfoToNormalize(fast_many[["info1"]], fast_many[["info2"]],
                                     normalize, treesIndependentInfo)

    ret <- .FloorNumericalNoise(treesIndependentInfo - msi - msi, treesIndependentInfo)
    ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                         infoInBoth = normInfo,
                         InfoInTree = SplitwiseInfo, Combine = "+")
    return(ret)
  }

  msi <- MatchingSplitInfo(tree1, tree2, normalize = FALSE, diag = FALSE,
                           reportMatching = reportMatching)

  treesIndependentInfo <- .MaxValue(tree1, tree2, SplitwiseInfo)
  normInfo <- if (!is.function(normalize)) treesIndependentInfo
  ret <- .FloorNumericalNoise(treesIndependentInfo - msi - msi, treesIndependentInfo)
  ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                       infoInBoth = normInfo,
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
