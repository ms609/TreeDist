#' @rdname TreeDistance
#' @export
MatchingSplitInfo <- function (tree1, tree2 = NULL, normalize = FALSE,
                               reportMatching = FALSE, diag = TRUE) {
  unnormalized <- CalculateTreeDistance(MatchingSplitInfoSplits, tree1,
                                        tree2, reportMatching)
  
  if (diag && is.null(tree2)) {
    unnormalized <- as.matrix(unnormalized)
    diag(unnormalized) <- SplitwiseInfo(tree1)
  }
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = SplitwiseInfo, Combine = .PairMean)
}

#' @rdname TreeDistance
#' @export
MatchingSplitInfoDistance <- function (tree1, tree2 = NULL, 
                                       normalize = FALSE,
                                       reportMatching = FALSE) {
  msi <- MatchingSplitInfo(tree1, tree2, normalize = FALSE, diag = FALSE,
                           reportMatching = reportMatching)
  
  treesIndependentInfo <- .MaxValue(tree1, tree2, SplitwiseInfo)
  ret <- treesIndependentInfo - msi - msi
  ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                       infoInBoth = treesIndependentInfo,
                       InfoInTree = SplitwiseInfo, Combine = '+')
  
  ret[ret < .Machine$double.eps^0.5] <- 0 # In case of floating point inaccuracy
  attributes(ret) <- attributes(msi)
  # Return:
  ret
}


#' @rdname TreeDistance
#' @export
MatchingSplitInfoSplits <- function (splits1, splits2,
                                     nTip = attr(splits1, 'nTip'),
                                     reportMatching = FALSE) {
  
  GeneralizedRF(splits1, splits2, nTip, cpp_msi_distance, maximize = TRUE,
                reportMatching = reportMatching)
}
