#' Approximate the Subtree Prune and Regraft distance
#'
#' `SPRDist()` approximates the \acronym{SPR} distance between
#' trees.
#' 
#' The function currently defaults to the heuristic method of
#' \insertCite{deOliveira2008;textual}{TreeDist}, which purports to provide an
#' upper bound on the \acronym{SPR} distance (though exceptions exist).
#' Other approximations 
#' \insertCite{@e.g. @Hickey2008, @Goloboff2008SPR, @Whidden2018}{TreeDist} are
#' not yet implemented.
#'
#' @template tree12ListParams
#' @param method Character specifying which method to use to approximate the
#' SPR distance.  Currently defaults to `"deOliveira"``.
#' `"Rogue"` implements an experimental method whose details are pending
#' publication; this function is under development, and may be modified or
#' removed without notice.  Once formally validated, it is anticipated that this
#' method will become the default.
#' @param symmetric Deprecated (redundant after fix of
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
SPRDist <- function(tree1, tree2 = NULL, method = "deOliveira", symmetric) {
  UseMethod("SPRDist")
}

#' @rdname SPRDist
#' @export
SPRDist.phylo <- function(tree1, tree2 = NULL, method = "deOliveira", symmetric) {
  if (is.null(tree2)) {
    NULL
  } else if (inherits(tree2, "phylo")) {
    .SPRFunc(method)(tree1, tree2)
  } else {
    vapply(tree2, .SPRFunc(method), double(1), tree1)
  }
}

.SPRFunc <- function(method) {
  switch(pmatch(tolower(gsub("\\s", "", method)),
                c("deoliveira2008", "rogue")),
         .SPRPairDeO, .SPRRogue)
}

#' @rdname SPRDist
#' @export
SPRDist.list <- function(tree1, tree2 = NULL, method = "deOliveira", symmetric) {
  if (is.null(tree2)) {
    PairwiseDistances(
      RootTree(RenumberTips(tree1, tree1), 1),
      .SPRFunc(method),
      check = FALSE
    )
  } else if (inherits(tree2, 'phylo')) {
    vapply(tree1, .SPRFunc(method), double(1), tree2)
  } else {
    vapply(tree2, SPRDist, double(length(tree1)), tree1, method = method)
  }
}

#' @rdname SPRDist
#' @export
SPRDist.multiPhylo <- SPRDist.list

# if x <- which(condition, arr.ind = FALSE), and
# y <-  which(condition, arr.ind = TRUE), then
# .Which1 and .Which2 give y[1] and y[2].
.Which1 <- function (x, nSplits) {
  ret <- x %% nSplits
  if (ret == 0L) {
    nSplits
  } else {
    ret
  }
}
.Which2 <- function (x, nSplits) (x - 1) %/% nSplits + 1L

.SPRExact5 <- function(sp1, sp2) {
  # Trees have shape (r, (p1, p2), (q1, q2)) after reduction.
  # There are two cases:
  if (all(xor(sp1[[1]] , sp1[[2]]) == xor(sp2[[1]] , sp2[[2]]))) {
    # Case 1: `r` has the same label in each tree
    # As the tree cannot be reduced by the reduction rule, we have
    # (X, (A, B), (C, D)) vs (X, (A, C), (B, C)) = 2 moves
    #
    2
  } else {
    # Case 2: `r` has a different label in each tree.
    # Assign `r` the label X in tree 1 and Y in tree 2.
    # (X, (Y, ?), (?, ?)) vs (Y, (X, ?), (?, ?))
    #
    # Then label the sister to X in tree2 X', and Y mutas mutantis
    # Notice that if X' == Y', the unlabelled cherry reduces by reduction rule
    # Hence we have
    # (X, ((Y, Y'), (X', ?))) vs (Y, ((X, X'), (Y', ?))) = 1 moves
    1
  }
}

.SPRExact6 <- function(sp1, sp2) {
  # Surprisingly, there is only one configuration with a distance of 1:
  # (((Lb, Lc), La), (Ra, (Rb, Rc))) vs (((Ra, Rb), La), (Lb, (Lc, Rc)))
  pairs1 <- TipsInSplits(sp1, smallest = TRUE) == 2
  pairs2 <- TipsInSplits(sp2, smallest = TRUE) == 2
  if (all(pairs1) || all(pairs2)) {
    return (2)
  }
  duo1 <- sp1[[pairs1]]
  trio1 <- sp1[[!pairs1]]
  
  .Overlapper <- function(s1, s2) {
    res <- as.logical(xor(s1, s2))
    if (sum(res) == 1) res else !res
  }
  middle1a <- .Overlapper(duo1[[1]], trio1)
  middle1b <- .Overlapper(duo1[[2]], trio1)
  
  duo2 <- sp2[[pairs2]]
  trio2 <- sp2[[!pairs2]]
  middle2a <- .Overlapper(duo2[[1]], trio2)
  middle2b <- .Overlapper(duo2[[2]], trio2)
  
  inMiddleEachTime <- (middle1a | middle1b) & (middle2a | middle2b)
  if (sum(inMiddleEachTime) == 1) {
    La <- if (inMiddleEachTime[middle1a]) middle1a else middle1b
    Lbc1 <- as.logical(duo1[[if (inMiddleEachTime[middle1a]) 1 else 2]])
    if (Lbc1[La]) {
      Lbc1 <- !Lbc1
    }
    if (any(Lbc1[middle2a | middle2b])) {
      # Lb is in the other middle position in tree 2:
      # (((?, ?), La), (Lb, (?, ?)))
      Lbc2 <- as.logical(duo2[[if (La[middle2a]) 1 else 2]])
      if (Lbc2[La]) {
        Lbc2 <- !Lbc2
      }
      if (!any(Lbc2[Lbc1])) {
        #   (((?, ?), La), (Lb, (Lc, ?)))
        #     (((Ra, Rb), La), (Lb, (Lc, Rc))) = 1 !!!
        return(1)
      }
      #   (((Lc, ?), La), (Lb, (?, ?)))
      #     (((Lc, Rc), La), (Lb, (Ra, Rb))) = 2
    }
    # (((?, ?), La), (Rb, (?, ?)))
    #   (((?, Ra), La), (Rb, (?, ?)))
    #     (((Lb, Ra), La), (Rb, (Rc, Lc))) = 2
    #   (((?, ?), La), (Rb, (Ra, ?)))
    #     (((Rc, Lc), La), (Rb, (Ra, Lb))) = 2
    # 
  }
  
  # All other tree pairs have a distance of 2 - see below
  return(2)
  
  # Trees may be one of two shapes: 
  #   ((a1, a2), (b1, b2), (c1, c2))
  #   (((Lb, Lc), La), (Ra, (Rb, Rc)))
  balanced1 <- all(pairs1)
  balanced2 <- all(pairs2)
  if (balanced1 && balanced2) {
    # There's only one possible configuration:
    # ((ab, ac), (ba, bc), (ca, cb)) vs ((ba, ca), (ab, cb), (ac, bc)) = 2
  }
  if (!balanced1 && !balanced2) {
    # Both trees have the shape
    # (((Lb, Lc), La), (Ra, (Rb, Rc)))
    # We will use the same labels for Tree 2, matching where possible.
    if (La1 == La2 && Ra1 == Ra2) {
      # La = La, Ra = Ra:
      # (((Lb, Lc), La), (Ra, (Rb, Rc))), (((Lb, Rb), La), (Ra, (Rc, Lc))) = 2
    }
    # As we can't match La and Ra, we'll match La if we can.
    if (La1 != La2 && Ra1 != Ra2) {
      # LO != La, Ra != Ra
      # La and Ra are both in the cherries
      # (((?, La), Lb), (Lc, (?, ?)))
      #   (((Rb, La), Lb), (Lc, (Ra, Rc))) = 2
      # 
      # (((?, ?), Rb), (Rc, (?, ?)))
      #   (((Lc, Ra), Rb), (Rc, (La, Lb))) = 2
      #   
      # (((?, ?), Lb), (Rb, (?, ?)))
      #   (((?, La), Lb), (Rb, (?, ?)))
      #     (((Lc, La), Lb), (Rb, (Ra, Rc))) = 2
      #     (((Ra, La), Lb), (Rb, (Lc, Rc))) = 2
      #     (((Rc, La), Lb), (Rb, (Ra, Lc))) = 2
      #     
      #   (((?, ?), Lb), (Rb, (La, ?)))
      #     (((Ra, Rc), Lb), (Rb, (La, Lc))) = 2
      #     (((Lc, Rc), Lb), (Rb, (La, Ra))) = 2
      #     (((Ra, Lc), Lb), (Rb, (La, Rc))) = 2
      #   
      #   (((?, ?), Lb), (Rb, (?, ?)))
      #     (((Lc, Rc), Lb), (Rb, (La, Ra))) = 2
      #     (((Lc, Ra), Lb), (Rb, (La, Rc))) = 2
      #     (((Ra, Rc), Lb), (Rb, (La, Lc))) = 2
      #     (((Ra, La), Lb), (Rb, (Rc, Lc))) = 2
      #     (((Rc, La), Lb), (Rb, (Ra, Lc))) = 2
      #     (((Lc, La), Lb), (Rb, (Ra, Rc))) = 2
    }
    # Else exactly one of the bridging leaves is the same; call this La.
    # 
    # (((?, ?), La), (Rb, (?, ?)))
    #   (((?, Ra), La), (Rb, (?, ?)))
    #     (((Lb, Ra), La), (Rb, (Rc, Lc))) = 2
    #   (((?, ?), La), (Rb, (Ra, ?)))
    #     (((Rc, Lc), La), (Rb, (Ra, Lb))) = 2
    # 
    # (((?, ?), La), (Lb, (?, ?)))
    #   (((Lc, ?), La), (Lb, (?, ?)))
    #     (((Lc, Rc), La), (Lb, (Ra, Rb))) = 2
    #   (((?, ?), La), (Lb, (Lc, ?)))
    #     (((Ra, Rb), La), (Lb, (Lc, Rc))) = 1 !!!
    # 
    # 
  }
}

.SPRExact7 <- function(sp1, sp2) {
  spr_table_7(sp1, sp2)
}

# Takes a 'Rogue' approach: finds the leaf that introduces the most conflict,
# and nixes it.
#' @importFrom TreeTools FirstMatchingSplit
.SPRRogue <- function(tree1, tree2, check = TRUE) {
  moves <- 0
  
  ProxyDistance <- switch(
    pmatch(toupper(getOption("sprProxy", "C")), c("C", "P", "Q", "R")),
    ClusteringInfoDist,
    PhylogeneticInfoDistance,
    function(x, y) Quartet::QuartetDivergence(Quartet::QuartetStatus(x, y)),
    RobinsonFoulds
  )
  
  reduced <- ReduceTrees(tree1, tree2, check = check)
  
  while (!is.null(reduced)) {
    
    tr1 <- reduced[[1]]
    tr2 <- reduced[[2]]
    nTip <- NTip(tr1)
    if (nTip == 4 && getOption("sprShortcut", Inf) >= 4) {
      return(moves + 1)
    }
    
    sp1 <- as.Splits(tr1)
    sp2 <- as.Splits(tr2, tr1)
    if (nTip == 5 && getOption("sprShortcut", Inf) >= 5) {
      return(moves + .SPRExact5(sp1, sp2))
    }
    if (nTip == 6 && getOption("sprShortcut", Inf) >= 6) {
      return(moves + .SPRExact6(sp1, sp2))
    }
    if (nTip == 7 && getOption("sprShortcut", Inf) >= 7) {
      return(moves + .SPRExact7(sp1, sp2))
    }
    
    firstMatchedSplit <- FirstMatchingSplit(sp1, sp2)
    if (!isFALSE(getOption("sprMatches")) && firstMatchedSplit > 0) {
      # At least one split exists in both trees
      subtips1 <- as.logical(sp1[[firstMatchedSplit]])
      subtips2 <- !subtips1
      
      # Add dummy tip as placeholder for other half of tree
      subtips1[!subtips1][[1]] <- TRUE

      # Repeat for other half-tree
      subtips2[!subtips2][[1]] <- TRUE

      moves1 <- .SPRRogue(
        KeepTipPostorder(tr1, subtips1),
        KeepTipPostorder(tr2, subtips1)
      )
      moves2 <- .SPRRogue(
        KeepTipPostorder(tr1, subtips2),
        KeepTipPostorder(tr2, subtips2)
      )
      return(moves + moves1 + moves2)
    }
    
    labels <- TipLabels(tr1)
    scores <- numeric(length(labels))
    blank <- rep_len(TRUE, length(labels))

    depth <- max(getOption("sprDepth", 1), 1)
    
    
    .ScoreWithout <- function(idx) {
      keep <- blank
      keep[idx] <- FALSE
      
      outcome <- keep_and_reduce(tr1, tr2, keep)
      if (is.null(outcome[[1]])) {
        return(-Inf)
      }
      
      oSpl <- as.Splits(outcome)
      firstMatch <- FirstMatchingSplit(oSpl[[1]], oSpl[[2]])
      
      if (firstMatch > 0) {
        subtips1 <- as.logical(oSpl[[1]][[firstMatch]])
        subtips2 <- !subtips1
        # Anchor to shared edge
        subtips1[!subtips1][[1]] <- TRUE
        subtips2[!subtips2][[1]] <- TRUE
        
        sub1 <- ReduceTrees(KeepTipPostorder(outcome[[1]], subtips1),
                            KeepTipPostorder(outcome[[2]], subtips1))
        
        sub2 <- ReduceTrees(KeepTipPostorder(outcome[[1]], subtips2),
                            KeepTipPostorder(outcome[[2]], subtips2))
        
        # Return:
        ProxyDistance(sub1[[1]], sub1[[2]]) +
          ProxyDistance(sub2[[1]], sub2[[2]])
        
      } else {
        # Return:
        ProxyDistance(outcome[[1]], outcome[[2]])
      }
    }
    for (i in seq_along(labels)) {
      scores[[i]] <- .ScoreWithout(i)
      if (!is.finite(scores[[i]])) break
    }
    if (any(!is.finite(scores))) {
      depth <- 1
    }
    if (depth > 1) {
      pairs <- combn(seq_along(labels), 2)
      nPairs <- dim(pairs)[[2]]
      pairScores <- double(nPairs)
      for (i in seq_len(nPairs)) {
        pairScores[[i]] <- .ScoreWithout(pairs[, i])
        if (!is.finite(pairScores[[i]])) break
      }
    }
    
    drop <- logical(length(labels))
    couldDrop <- scores == min(scores)

    if (depth > 1 && min(pairScores) < min(scores)) {
      pairDrop <- pairScores == min(pairScores)
      dropTions <- pairs[, pairDrop]
      if (any(couldDrop[dropTions]) && !any(!is.finite(pairScores))) {
        # Dropping two at once doesn't give us any benefit over dropping one
        # at a time – but will mean we can't spot a handy pair next time.
        drop[dropTions[which.min(scores[dropTions])]] <- TRUE
      } else {
        # If dropping two gives us a better solution, drop both at once -
        # failing to do so can cause an optimal path to be missed.
        drop[pairs[, which.max(pairDrop)]] <- TRUE
      }
    } else {
      drop[[which.min(scores)]] <- TRUE
    }
    
    reduced <- keep_and_reduce(tr1, tr2, !drop)
    if (length(reduced) == 1L) {
      reduced <- NULL
    }
    
    moves <- sum(moves, drop)
  }
  
  # Return:
  moves
}

# An attempt to reproduce the phangorn results using the algorithm of 
# \insertCite{deOliveira2008;textual}{TreeDist}
# An exact match is unlikely due to the arbitrary breaking of ties when leaves
# are selected for removal.
#' @examples
#' # de Oliveira Martins et al 2008, fig. 7
#' tree1 <- ape::read.tree(text = "((1, 2), ((a, b), (c, d)), (3, (4, (5, (6, 7)))));")
#' tree2 <- ape::read.tree(text = "((1, 2), 3, (4, (5, (((a, b), (c, d)), (6, 7)))));")
#' oPar <- par(mfrow =c(2, 1), mar = rep(0, 4))
#' plot(tree1)
#' plot(tree2)
#' par(oPar)
#' SPRDist(tree1, tree2, method = "deO")
#' @keywords internal
#' @importFrom TreeTools DropTip TipsInSplits KeepTipPostorder
#' @importFrom TreeTools edge_to_splits
.SPRPairDeO <- function(tree1, tree2, check = TRUE) {
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
      unmatchedSplits <- is.na(matched[["matching"]])
      sp1 <- sp1[[unmatchedSplits]]
      sp2 <- sp2[[-matched$matching[!unmatchedSplits]]]
    }
    
    nSplits <- length(sp1)
    # Compute size of disagreement splits - see Fig. 7C in @deOliv2008
    mmSize <- mismatch_size(sp1, sp2)
    stopifnot(all(mmSize > 0))
    
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
