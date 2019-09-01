#' @describeIn TreeDistance Mutual Matching Split information of two trees.
#' @export
MutualMatchingSplitInfo <- function (tree1, tree2, normalize = FALSE, 
                                     reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(MutualMatchingSplitInfoSplits, tree1, tree2,
                                        reportMatching)
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = PartitionInfo, Combine = pmin)
}

#' @describeIn TreeDistance Variation of matching split information between two trees.
#' @export
VariationOfMatchingSplitInfo <- function (tree1, tree2, normalize = FALSE,
                                          reportMatching = FALSE) {
  mmsi <- MutualMatchingSplitInfo(tree1, tree2, normalize = FALSE, reportMatching)
  
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


#' @describeIn TreeDistance Calculate variation of matching split information from splits instead of trees.
#' @inheritParams MutualPhylogeneticInfoSplits
#' @importFrom TreeSearch LnUnrooted.int LogTreesMatchingSplit
#' @export
MutualMatchingSplitInfoSplits <- function (splits1, splits2, reportMatching = FALSE) {
  dimSplits1 <- dim(splits1)
  dimSplits2 <- dim(splits2)
  nTerminals <- dimSplits1[1]
  if (dimSplits2[1] != nTerminals) {
    stop("Split rows must bear identical labels")
  }
  
  swapSplits <- (dimSplits1[2] > dimSplits2[2])
  if (swapSplits) {
    # solve_LDAP expects splits1 to be no larger than splits2
    tmp <- splits1
    splits1 <- splits2
    splits2 <- tmp
    
    tmp <- dimSplits1
    dimSplits1 <- dimSplits2
    dimSplits2 <- tmp
    
    remove(tmp)
  }
  
  taxonNames1 <- rownames(splits1)
  taxonNames2 <- rownames(splits2)
  
  if (!is.null(taxonNames2)) {
    splits2 <- unname(splits2[taxonNames1, , drop=FALSE])
    splits1 <- unname(splits1) # split2[split1] faster without names
  }
  
  nSplits1 <- dimSplits1[2]
  nSplits2 <- dimSplits2[2]
  if (nSplits1 == 0) return (0)
  
  AgreementInfoNats <- function (splitI, agree) {
    inAgreement <- sum(agree)
    n0 <- sum(splitI[agree])
    LnUnrooted.int(inAgreement) - LogTreesMatchingSplit(n0, inAgreement - n0)
  }
  
  pairScores <- matrix((mapply(function(i, j) {
    splitI0 <- splits1[, i]
    splitJ0 <- splits2[, j]
    
    agree1 <- splitI0 == splitJ0
    
    max(AgreementInfoNats(splitI0, agree1),
        AgreementInfoNats(splitI0, !agree1))
    
  }, seq_len(nSplits1), rep(seq_len(nSplits2), each=nSplits1)
  )), nSplits1, nSplits2) / log(2)
  
  TreeDistanceReturn(pairScores, maximize = TRUE,
                     reportMatching, swapSplits,
                     splits1, splits2,
                     taxonNames1)
}
