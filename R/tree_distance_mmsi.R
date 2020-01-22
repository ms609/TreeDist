#' @describeIn TreeDistance Mutual Matching Split information of two trees.
#' @export
MatchingSplitInfo <- function (tree1, tree2 = tree1, normalize = FALSE,
                                     reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(MatchingSplitInfoSplits, tree1,
                                        tree2, reportMatching)
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = SplitwiseInfo, Combine = pmin)
}

#' @describeIn TreeDistance Variation of matching split information between two
#' trees.
#' @export
MatchingSplitInfoDistance <- function (tree1, tree2 = tree1, 
                                          normalize = FALSE,
                                          reportMatching = FALSE) {
  msi <- MatchingSplitInfo(tree1, tree2, normalize = FALSE, 
                                  reportMatching)
  
  treesIndependentInfo <- outer(SplitwiseInfo(tree1), SplitwiseInfo(tree2), '+')
  ret <- treesIndependentInfo - msi - msi
  ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                       infoInBoth = treesIndependentInfo,
                       InfoInTree = SplitwiseInfo, Combine = '+')
  
  ret[ret < 1e-13] <- 0 # In case of floating point inaccuracy
  attributes(ret) <- attributes(msi)
  # Return:
  ret
}


#' @rdname TreeDistance
#' @export
MatchingSplitInfoSplits <- function (splits1, splits2,
                                           nTip = attr(splits1, 'nTip'),
                                           reportMatching = FALSE) {
  
  GeneralizedRF(splits1, splits2, nTip, cpp_mmsi_distance, maximize = TRUE,
                reportMatching = reportMatching)
}
