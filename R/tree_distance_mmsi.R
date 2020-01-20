#' @describeIn TreeDistance Mutual Matching Split information of two trees.
#' @export
MutualMatchingSplitInfo <- function (tree1, tree2 = tree1, normalize = FALSE,
                                     reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(MutualMatchingSplitInfoSplits, tree1,
                                        tree2, reportMatching)
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = PartitionInfo, Combine = pmin)
}

#' @describeIn TreeDistance Variation of matching split information between two
#' trees.
#' @export
VariationOfMatchingSplitInfo <- function (tree1, tree2 = tree1, 
                                          normalize = FALSE,
                                          reportMatching = FALSE) {
  mmsi <- MutualMatchingSplitInfo(tree1, tree2, normalize = FALSE, 
                                  reportMatching)
  
  treesIndependentInfo <- outer(PartitionInfo(tree1), PartitionInfo(tree2), '+')
  ret <- treesIndependentInfo - mmsi - mmsi
  ret <- NormalizeInfo(ret, tree1, tree2, how = normalize,
                       infoInBoth = treesIndependentInfo,
                       InfoInTree = PartitionInfo, Combine = '+')
  
  ret[ret < 1e-13] <- 0 # In case of floating point inaccuracy
  attributes(ret) <- attributes(mmsi)
  # Return:
  ret
}


#' @describeIn TreeDistance Calculate variation of matching split information
#'   from splits instead of trees.
#' @inheritParams SharedPhylogeneticInfoSplits
#' @export
MutualMatchingSplitInfoSplits <- function (splits1, splits2,
                                           nTip = attr(splits1, 'nTip'),
                                           reportMatching = FALSE) {
  
  GeneralizedRF(splits1, splits2, nTip, cpp_mmsi_distance, maximize = TRUE,
                reportMatching = reportMatching)
}
