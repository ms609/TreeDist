#' Approximate the Subtree Prune and Regraft (SPR) distance.
#' 
#' `SPRDist()` calculates an upper bound on the SPR distance between trees
#' using the heuristic method of \insertCite{deOliveira2008;textual}{TreeDist}.
#' Other approximations are available
#' \insertCite{@e.g. @Hickey2008, @Goloboff2008SPR, @Whidden2018}{TreeDist}.
#' 
#' @template tree12ListParams
#' @param method Character specifying which method to use to approximate the
#' SPR distance.  Currently defaults to `"deOliveira"`, the only available
#' option; a new method will eventually become the default.
#' @param symmetric Ignored (redundant after fix of
#' [phangorn#97](https://github.com/KlausVigo/phangorn/issues/97)).
#' 
#' @return `SPRDist()` returns a vector or distance matrix of distances 
#' between trees.
#' 
#' @references \insertAllCited{}
#' 
#' @examples
#' library("TreeTools", quietly = TRUE)
#' 
#' # Compare single pair of trees
#' SPRDist(BalancedTree(7), PectinateTree(7))
#' 
#' # Compare all pairs of trees        
#' SPRDist(as.phylo(30:33, 8))
#' 
#' # Compare each tree in one list with each tree in another
#' SPRDist(BalancedTree(7), as.phylo(0:2, 7))
#' SPRDist(as.phylo(0:2, 7), PectinateTree(7))
#'
#' SPRDist(list(bal = BalancedTree(7), pec = PectinateTree(7)),
#'         as.phylo(0:2, 7))
#' @template MRS
#'   
#' @seealso Exact calculation with [\pkg{TBRDist}](
#' https://ms609.github.io/TBRDist/reference/TreeRearrangementDistances.html)
#' functions `USPRDist()` and `ReplugDist()`.
#' 
#' \pkg{phangorn} function \code{\link[phangorn:treedist]{SPR.dist()}} employs
#' the \insertCite{deOliveira2008;textual}{TreeDist} algorithm but can crash
#' when sent trees of certain formats, and tends to have a longer running time.
#' 
#' @family tree distances
#' @importFrom TreeTools PairwiseDistances Postorder
#' @export
SPRDist <- function (tree1, tree2 = NULL, method = "deOliveira", symmetric) {
  UseMethod("SPRDist")
}

#' @rdname SPRDist
#' @export
SPRDist.phylo <- function (tree1, tree2 = NULL, method = "deOliveira", symmetric) {
  if (is.null(tree2)) {
    NULL
  } else if (inherits(tree2, "phylo")) {
    .SPRFunc(method)(tree1, tree2)
  } else {
    vapply(tree2, .SPRFunc(method), double(1), tree1)
  }
}

.SPRFunc <- function(method) {
  .SPRPairDeO
}

#' @rdname SPRDist
#' @export
SPRDist.list <- function (tree1, tree2 = NULL, method = "deOliveira", symmetric) {
  if (is.null(tree2)) {
    PairwiseDistances(RootTree(RenumberTips(tree1, tree1), 1),
                      .SPRFunc(method), check = FALSE)
  } else if (inherits(tree2, 'phylo')) {
    vapply(tree1, .SPRFunc(method), double(1), tree2)
  } else {
    vapply(tree2, SPRDist, double(length(tree1)), tree1, method = method)
  }
}

#' @rdname SPRDist
#' @export
SPRDist.multiPhylo <- SPRDist.list

.Which1 <- function (x, nSplits) {
  ret <- x %% nSplits
  if (ret == 0L) {
    nSplits
  } else {
    ret
  }
}
.Which2 <- function (x, nSplits) (x - 1) %/% nSplits + 1L

#' @importFrom TreeTools DropTip TipsInSplits KeepTipPostorder
#' @importFrom TreeTools edge_to_splits
.SPRPairDeO <- function(tree1, tree2, check = TRUE, debug = FALSE) {
  moves <- 0
  
  # Reduce trees (Fig. 7A in deOliveira2008)
  reduced <- ReduceTrees(tree1, tree2, check = check)
  
  while (!is.null(reduced)) {
    tr1 <- reduced[[1]]
    tr2 <- reduced[[2]]
    edge1 <- tr1[["edge"]]
    edge2 <- tr2[["edge"]]
    labels <- tr1[["tip.label"]]
    nTip <- length(labels)
    sp1 <- edge_to_splits(edge1, PostorderOrder(edge1), labels, nTip = nTip)
    sp2 <- edge_to_splits(edge2, PostorderOrder(edge2), labels, nTip = nTip)
    matched <- cpp_robinson_foulds_distance(sp1, sp2, nTip)
    nMatched <- matched[["score"]]
    if (nMatched != length(sp1) * 2) {
      if (debug) {
        message("Identical splits: ", length(sp1) - (nMatched / 2))
      }
      unmatchedSplits <- is.na(matched$matching)
      sp1 <- sp1[[unmatchedSplits]]
      sp2 <- sp2[[-matched$matching[!unmatchedSplits]]]
    }
    
    nSplits <- length(sp1)
    # Compute size of disagreement splits - see Fig. 7C in @deOliv2008
    mmSize <- mismatch_size(sp1, sp2)
    if (any(mmSize == 0)) {
      message("Zero-sizers!")
      sapply(which(mmSize == 0), .Which1, nSplits)
      sapply(which(mmSize == 0), .Which2, nSplits)
    }
    
    # Arbitrary selection of leaves to remove introduces a stochastic element
    minMismatch <- which.min(mmSize)
    
    split1 <- structure(sp1[.Which1(minMismatch, nSplits), , drop = FALSE],
                        nTip = nTip, class = "Splits")
    split2 <- structure(sp2[.Which2(minMismatch, nSplits), , drop = FALSE],
                        nTip = nTip, class = "Splits")
    disagreementSplit <- xor(split1, split2)
    
    keep <- as.logical(disagreementSplit)
    nKeep <- sum(keep)
    if (nKeep < length(keep) / 2) {
      keep <- !keep
      nKeep <- length(keep) - nKeep
    }
    
    reduced <- if (nKeep < 4L) {
      NULL
    } else {
      keep_and_reduce(tr1, tr2, keep)
    }
    
    if (length(reduced) == 1L) {
      reduced <- NULL
    }
    
    moves <- moves + 1 # Usually an underestimate, unless we've lost a chance
                       # to "untangle a knot"
  }
  
  # Return:
  moves
}
