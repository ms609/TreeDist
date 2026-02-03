#' Approximate the Subtree Prune and Regraft distance
#'
#' `SPRDist()` calculates an upper bound on the SPR distance between trees
#' using the heuristic method of \insertCite{deOliveira2008;textual}{TreeDist}.
#' Other approximations are available
#' \insertCite{@e.g. @Hickey2008, @Goloboff2008SPR, @Whidden2018}{TreeDist}.
#'
#' @template tree12ListParams
#' @param method Character specifying which method to use to approximate the
#' SPR distance.  Currently defaults to `"deOliveira"``, the only available
#' option; a new method will become the default once available.
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
SPRDist <- function(tree1, tree2 = NULL, method = "confl", symmetric) {
  UseMethod("SPRDist")
}

#' @rdname SPRDist
#' @export
SPRDist.phylo <- function(tree1, tree2 = NULL, method = "confl", symmetric) {
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
                c("deoliveira2008", "confl", "experiment")),
         .SPRPairDeOCutter, .SPRConfl, .SPRExperiment)
}

#' @rdname SPRDist
#' @export
SPRDist.list <- function(tree1, tree2 = NULL, method = "confl", symmetric) {
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

.SPRExperiment <- function(tree1, tree2, check = TRUE) {
  debug <- isTRUE(getOption("debugSPR", FALSE))
  moves <- 0
  if (debug) dropList <- character(0)
  
  reduced <- ReduceTrees(tree1, tree2, check = check)
  if (debug) {
    dropList <- character(0)
    par(mfrow = 1:2, mai = rep(0.1, 4))
    oldBG <- par(bg = "#eeddcc")
    plot(reduced[[1]])
    nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
    plot(reduced[[2]])
    nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
    par(oldBG)
  }
  
  while (!is.null(reduced)) {
    tr1 <- reduced[[1]]
    tr2 <- reduced[[2]]
    edge1 <- tr1[["edge"]]
    edge2 <- tr2[["edge"]]
    labels <- tr1[["tip.label"]]
    nTip <- length(labels)
    sp1 <- edge_to_splits(edge1, PostorderOrder(edge1), labels, nTip = nTip)
    sp2 <- edge_to_splits(edge2, PostorderOrder(edge2), labels, nTip = nTip)
    nSplits <- length(sp1)
    
    confusion <- confusion(sp1, sp2)
    if (debug) {
      dimnames(confusion) <- list(
        c("ab", "aB", "Ab", "AB"),
        names(sp1),
        names(sp2))
    }
    concave <- colSums(confusion == 0)
    
    matches <- concave == 2
    if (any(matches)) {
      agreement <- as.logical(sp1[[.Which1(which.max(matches), nSplits)]])
      if (debug) {
        action <- paste("Two identical subtrees",
                        names(sp1[[.Which1(which.max(matches), nSplits)]]),
                        " = ",
                        names(sp2[[.Which2(which.max(matches), nSplits)]])
        )
        nodelabels("|______", frame = "n", col = "darkred", font = 3,
                   as.integer(names(sp2[[.Which2(which.max(matches), nSplits)]])))
        legend("topleft", action, bty = "n")
        message(action)
      }
      subtips1 <- agreement
      subtips1[!subtips1][1] <- TRUE
      subtips2 <- !agreement
      subtips2[agreement][1] <- TRUE
      return(
        moves +
          .SPRExperiment(
            KeepTipPostorder(tr1, subtips1),
            KeepTipPostorder(tr2, subtips1)
          ) +
          .SPRExperiment(
            KeepTipPostorder(tr1, subtips2),
            KeepTipPostorder(tr2, subtips2)
          )
      )
    }
    .Is1 <- function (i, j) {
      hitHere <- logical(attr(sp1, "nTip"))
      if (confusion[1, i, j] == 1) {
        hitHere <- hitHere | as.logical(sp1[[i]] & sp2[[j]])
      }
      if (confusion[2, i, j] == 1) {
        hitHere <- hitHere | as.logical(sp1[[i]] & !sp2[[j]])
      }
      if (confusion[3, i, j] == 1) {
        hitHere <- hitHere | as.logical(!sp1[[i]] & sp2[[j]])
      }
      if (confusion[4, i, j] == 1) {
        hitHere <- hitHere | as.logical(!sp1[[i]] & !sp2[[j]])
      }
      hitHere
    }
    
    .FindDrops <- function (x) {
      which(.Is1(.Which1(x, nSplits), .Which2(x, nSplits)))
    }
    
    .FindOverlap <- function (x) {
      i <- .Which1(x, nSplits)
      j <- .Which2(x, nSplits)
      which(
      if (confusion[1, i, j] > 1) {
        as.logical(!sp1[[i]] & !sp2[[j]])
      } else 
      if (confusion[2, i, j] > 1) {
        as.logical(!sp1[[i]] & sp2[[j]])
      } else
      if (confusion[3, i, j] > 1) {
        as.logical(sp1[[i]] & !sp2[[j]])
      } else
      if (confusion[4, i, j] > 1) {
        as.logical(sp1[[i]] & sp2[[j]])
      })
    }
    
    nits <- which(apply(confusion, 2:3, function (x) sum(0:2 %in% x)) == 3)
    nitDrops <- vapply(nits, function (x) which(.Is1(.Which1(x, nSplits), .Which2(x, nSplits))), integer(1))
    nitDups <- duplicated(nitDrops)
    if (any(nitDups)) {
      if (debug) {
        message("Doubly dropped: ", paste(nitDrops[nitDups], collapse = ", "))
      }
      nitDrops <- nitDrops[nitDups]
      #nitDrops <- nitDrops[!nitDups]
    }
    if (debug) {
      if (length(nits)) {
        i <- sapply(nits, .Which1, nSplits)
        j <- sapply(nits, .Which2, nSplits)
        message(paste(
          apply(cbind(names(sp1)[i], names(sp2)[j]), 1, paste, collapse = "-"),
          collapse = "; ")
        )
        message(" -> Drop options: ", paste(labels[nitDrops], collapse = ", "))
      } else {
        message(" -> No nit options")
      }
    }
    twits <- double(0)
    if (!length(nitDrops)) {
      twits <- which(apply(confusion, 2:3, function (x) sum(x == 1) > 1))
      if (length(twits)) {
        twitDrops <- unlist(sapply(twits, .FindDrops))
        keep <- !tabulate(which.max(tabulate(twitDrops)), nTip)
      
        # flits <- which(apply(confusion, 2:3, function (x) sum(x == 1) == 3))
        # flitDrops <- vapply(flits, .FindOverlap, integer(1))
        # nitDrops <- unique(flitDrops)
        # if (debug) {
        #   message("    Flit candidates: ",
        #           paste(labels[nitDrops], collapse = ", "),
        #           "; duplicates: ",
        #           paste(labels[flitDrops[duplicated(as.integer(flitDrops))]],
        #                 collapse = ", ")
        #   )
        # }
      }
    }
    if (!length(twits)) {
      hits <- double(length(labels))
      for (i in seq_along(sp1)) for (j in seq_along(sp2)) {
        hitHere <- .Is1(i, j)
        hits[hitHere] <- hits[hitHere] + 1L
      }
      if (length(nitDrops)) {
        nitHits <- hits[nitDrops]
        keep <- !tabulate(min(nitDrops[nitHits == max(nitHits)]), nTip)
        
        if (debug) {
          message(" -> Most-hit nit = ",
                  labels[nitDrops[nitHits == max(nitHits)]],
                  " (", max(nitHits),"); cf. ",
                  labels[hits == max(hits)], " (", max(hits), ")")
        }
      } else {
        if (debug) {
          message(" -> Sore thumbs: ",
                  paste(labels[hits == max(hits)], collapse = ", "))
        }
      
        keep <- !tabulate(which.max(hits), length(labels))
      }
    }
    nKeep <- sum(keep)
    if (debug) {
      drop <- !keep
      dropList <- c(dropList, labels[drop])
      action <- paste0("Dropping: ", paste(labels[drop], collapse = ", "),
                       " (", paste(which(drop), collapse = ", "), ")")
      legend("topleft", action, bty = "n")
      ape::tiplabels("|________", which(drop),
                     frame = "n", col = "darkred", font = 2)
      message(action)
    }
    reduced <- if (nKeep < 4L) {
      NULL
    } else {
      keep_and_reduce(tr1, tr2, keep)
    }
    if (length(reduced) == 1) {
      reduced <- NULL
    }
    
    if (debug) {
      if (is.null(reduced)) {
        plot.new(); plot.new()
      } else {
        plot(reduced[[1]])
        nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
        plot(reduced[[2]])
        nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
      }
    }
    
    moves <- moves + 1
  }
  
  # Return:
  
  if (debug) list(moves, dropList) else moves
}

# tree1, tree2 are phylo objects
# ReduceTrees implements the reduction rule, i.e. collapsing identical cherries
# confusion() gives the confusion matrix of tips in each partition of two splits
.SPRConfl <- function(tree1, tree2, check = TRUE) {
  moves <- 0
  debug <- isTRUE(getOption("debugSPR", FALSE))
  if (debug) dropList <- character(0)
  
  reduced <- ReduceTrees(tree1, tree2, check = check)
  if (!is.null(reduced) && debug) {
    par(mfrow = 1:2, mai = rep(0.1, 4))
    plot(reduced[[1]])
    ape::nodelabels(frame = "none", cex = 0.8)
    plot(reduced[[2]])
    ape::nodelabels(frame = "none", cex = 0.8)
  }
  
  while (!is.null(reduced)) {
    tr1 <- reduced[[1]]
    tr2 <- reduced[[2]]
    edge1 <- tr1[["edge"]]
    edge2 <- tr2[["edge"]]
    labels <- tr1[["tip.label"]]
    nTip <- length(labels)
    sp1 <- edge_to_splits(edge1, PostorderOrder(edge1), labels, nTip = nTip)
    sp2 <- edge_to_splits(edge2, PostorderOrder(edge2), labels, nTip = nTip)
    nSplits <- length(sp1)
    
    confusion <- confusion(sp1, sp2)
    
    concave <- colSums(confusion == 0)
    
    # Divide and conquer can help - but doesn't always.
    # 
    # If we have an edge that's shared in both trees, we know that's now
    # part of our best possible solution . The question becomes: where do we
    # want that edge?  
    # If we relocate that edge smartly, we can greatly increase the number
    # of common edges.  "Knots" are untangled when we increase the number of 
    # common edges by lots.
    # 
    # Consider the pectinate knot:
    #    z-y-x-...-h-g=F=D-E-[cba]
    #    h-g-i-...-y-z=F=E-D-[cba]
    #    
    # We have two edges in common (=).  If we break off |f d e cba, and attach
    # next to h, then we gain OODLES of common edges and can massively reduce.
    # 
    # If we cut at g=F / z=F, then we have one move to detach the anchor in
    # z-y-x-...-h-g-*
    # h-g-i-...-y-z-*
    # 
    # Then in our other half we can excise any of D, E, cba to reconcile
    # [*F]-D-E-cba
    # [*F]-E-D-cba
    # 
    # Without this subdivision we have to lop our way through cba, E, D to reach
    # the common edge.
    # 
    matches <- concave == 2
    if (!isFALSE(getOption("sprMatches")) && any(matches)) {
      # At least one split exists in both trees
      matchingSplit <- which.max(matches)
      agreement <- as.logical(sp1[[.Which1(matchingSplit, nSplits)]])

      # Take left side of split
      subtips1 <- agreement
      # Add dummy tip as placeholder for other half of tree
      subtips1[!subtips1][[1]] <- TRUE

      # Repeat for other half-tree
      subtips2 <- !agreement
      subtips2[!subtips2][[1]] <- TRUE

      if (debug) {
        message("Division A: ", 
                paste(colnames(agreement)[agreement], collapse = " "), 
                " | ",
                paste(colnames(agreement)[!agreement], collapse = " "))
        colNow <- par("col")
        if (colNow == "black") colNow <- "#000000"
        colIdx <- match(colNow, palette.colors(8), 0)
        oPar <- par(col = palette.colors(8)[colIdx + 1])
        on.exit(par(oPar))
      }
      moves1 <- .SPRConfl(
        KeepTipPostorder(tr1, subtips1),
        KeepTipPostorder(tr2, subtips1)
      )
      if (debug) {
        message("Division B: ", paste(colnames(agreement)[!agreement], collapse = " "))
        colNow <- par("col")
        colIdx <- match(colNow, palette.colors(8), 0)
        par(col = palette.colors(8)[colIdx + 1])
      }
      moves2 <- .SPRConfl(
        KeepTipPostorder(tr1, subtips2),
        KeepTipPostorder(tr2, subtips2)
      )
      return(if (debug) {
        structure(
          moves + moves1 + moves2,
          dropList = paste(
            dropList,
            attr(moves1, "dropList"),
            attr(moves2, "dropList"),
            collapse = " | ", sep = " ")
          )
        } else {
        moves + moves1 + moves2
      })
    }
    
    confInf <- confusion
    confInf[confusion == 0] <- Inf
    confMin <- apply(confInf, 2:3, min)
    minConf <- min(confMin[confMin > 0])
    if (debug && minConf > 1) {
      message("Minimum conflict: ", minConf)
    }
    minConfOpts <- confMin == minConf
    
    candidates <- switch(
      pmatch(tolower(getOption("sprH", "confusion")), c("confusion", "vi",
                                                        "ami", "joint", "vinorm")),
      {
        
        # if minConf == 1, then removing a single leaf can resolve a contradiction
        # We use entropy to decide which leaf might be most profitable to remove.
        # 
        # h gives the joint entropy of each pair of splits in tree1 & tree2
        h <- apply(confusion, 2:3, Ntropy)
        minH <- min(h[minConfOpts])
        maxH <- max(h[minConfOpts])
        
        candidates <- which(h == maxH, arr.ind = TRUE)
        if (!isFALSE(getOption("sprTies")) && nrow(candidates) > 1) {
          # Let's identify the split in each tree that is most at odds with all
          # other splits
          tieBreak <- outer(rowMeans(h), colMeans(h))[candidates]
          candidates <- candidates[tieBreak == max(tieBreak), , drop = FALSE]
        }
        
        # If still tied, break arbitrarily.
        if (nrow(candidates) > 1) {
          # TODO perhaps we can find a non-arbitrary way to break any remaining ties?
          neyms <- cbind(names(sp1)[candidates[, 1]], names(sp2)[candidates[, 2]])
          message("Candidates remain tied: ",
                  paste(apply(neyms, 1, paste, collapse = "-"), collapse = ", "))
        }
        candidates
      },
      { # vi
        # In AZ33, gets "distracted" by Z, Y, L before finding the AD trick
        nTip <- NTip(sp1)
        in1 <- TipsInSplits(sp1)
        in2 <- TipsInSplits(sp2)
        n1 <- rbind(in1, nTip - in1)
        n2 <- rbind(in2, nTip - in2)
        
        h1 <- apply(n1, 2, Ntropy)
        h2 <- apply(n2, 2, Ntropy)
        h12 <- apply(confusion, 2:3, Ntropy)
        mi <- outer(h1, h2, "+") - h12
        
        emi <- outer(seq_along(in1), seq_along(in2),
                     Vectorize(function(i, j) expected_mi(n1[, i], n2[, j])))
        ami <- mi - emi
        vi <- outer(h1, h2, "+") - (ami + ami)
        vi[!minConfOpts] <- -Inf
        which(vi == max(vi), arr.ind = TRUE)
      },
      { # ami
        nTip <- NTip(sp1)
        in1 <- TipsInSplits(sp1)
        in2 <- TipsInSplits(sp2)
        n1 <- rbind(in1, nTip - in1)
        n2 <- rbind(in2, nTip - in2)
        
        h1 <- apply(n1, 2, Ntropy)
        h2 <- apply(n2, 2, Ntropy)
        h12 <- apply(confusion, 2:3, Ntropy)
        mi <- outer(h1, h2, "+") - h12
        
        emi <- outer(seq_along(in1), seq_along(in2),
                     Vectorize(function(i, j) expected_mi(n1[, i], n2[, j])))
        ami <- mi - emi
        score <- ami
        score[!minConfOpts] <- Inf
        which(score == min(score), arr.ind = TRUE)
      },
      { # joint
        nTip <- NTip(sp1)
        in1 <- TipsInSplits(sp1)
        in2 <- TipsInSplits(sp2)
        n1 <- rbind(in1, nTip - in1)
        n2 <- rbind(in2, nTip - in2)
        
        h1 <- apply(n1, 2, Ntropy)
        h2 <- apply(n2, 2, Ntropy)
        h12 <- apply(confusion, 2:3, Ntropy)
        score <- h12
        score[!minConfOpts] <- -Inf
        which(score == max(score), arr.ind = TRUE)
      },
      { # viNorm
        nTip <- NTip(sp1)
        in1 <- TipsInSplits(sp1)
        in2 <- TipsInSplits(sp2)
        n1 <- rbind(in1, nTip - in1)
        n2 <- rbind(in2, nTip - in2)
        
        h1 <- apply(n1, 2, Ntropy)
        h2 <- apply(n2, 2, Ntropy)
        h12 <- apply(confusion, 2:3, Ntropy)
        mi <- outer(h1, h2, "+") - h12
        
        emi <- outer(seq_along(in1), seq_along(in2),
                     Vectorize(function(i, j) expected_mi(n1[, i], n2[, j])))
        ami <- mi - emi
        vi <- outer(h1, h2, "+") - (ami + ami)
        
        score <- vi / outer(h1, h2, "+")
        candidates <- which(score == max(score[minConfOpts]) & minConfOpts,
                            arr.ind = TRUE)
        if (nrow(candidates) > 1) {
          tieBreak <- outer(rowMeans(score), colMeans(score))[candidates]
          candidates <- candidates[tieBreak == max(tieBreak), , drop = FALSE]
        }
      }
    )
    
    splitA <- sp1[[candidates[1, 1]]]
    splitB <- sp2[[candidates[1, 2]]]
    ins <- TipsInSplits(c(splitA, splitB, splitA & splitB),
                        keep.names = FALSE)
    nTip <- attr(splitA, "nTip")
    AB <- ins[3]
    aB <- ins[2] - ins[3]
    Ab <- ins[1] - ins[3]
    ab <- nTip - (aB + Ab + AB)
    confusion <- c(AB = AB, ab = ab, aB = aB, Ab = Ab)
    if (debug) {
      summary(c(splitA, splitB))
      print(confusion)
    }
    balance <- ins[1:2] + ins[1:2] - nTip
    absBal <- abs(balance)
    # balance > 0 means more in than out
    confusInf <- ifelse(confusion == 0, Inf, confusion)
    mins <- confusion == min(confusInf)
    drop <- switch(sum(mins),
      switch(which(mins),
             as.logical(splitA & splitB),
             as.logical(!splitA & !splitB),
             as.logical(!splitA & splitB),
             as.logical(splitA & !splitB)),
      {
        if (absBal[1] > absBal[2]) {
          # remove from A, which is less balanced
          if (balance[1] < 0) {
            # remove from IN A
            if (confusion[1] == confusion[4]) {
              as.logical(splitA)
            } else if (confusInf[1] < confusInf[4]) {
              as.logical(splitA & splitB)
            } else {
              as.logical(splitA & !splitB)
            }
          } else {
            # remove from OUT A
            if (confusion[2] == confusion[3]) {
              as.logical(!splitA)
            } else if (confusInf[2] < confusInf[3]) {
              as.logical(!splitA & !splitB)
            } else {
              as.logical(!splitA & splitB)
            }
          }
        } else {
          # remove from B, which is less balanced
          if (balance[2] < 0) {
            # remove from IN B
            if (confusion[1] == confusion[3]) {
              as.logical(splitB)
            } else if (confusInf[1] < confusInf[3]) {
              as.logical(splitA & splitB)
            } else {
              as.logical(!splitA & splitB)
            }
          } else {
            # remove from OUT B
            if (confusion[2] == confusion[4]) {
              as.logical(!splitB)
            } else if (confusInf[2] < confusInf[4]) {
              as.logical(!splitA & !splitB)
            } else {
              as.logical(splitA & !splitB)
            }
          }
        }
      }, 
      if (abs(Ab - aB) > abs(AB - ab)) {
        if (aB < Ab) {
          as.logical(!splitA & splitB)
        } else {
          as.logical(splitA & !splitB)
        }
      } else {
        if (ab < AB) {
          as.logical(!splitA & !splitB)
        } else {
          as.logical(splitA & splitB)
        } 
      }, 
      if (sum(mins) == 4) {
        1:4 == 1
      } else {
        print(confusion)
        summary(c(splitA, splitB))
        stop()
      })

    if (!any(drop)) {
      splitA <<- splitA
      splitB <<- splitB
      stopifnot(any(drop))
    }
    if (sum(drop) > 1) {
      drop <- as.logical(tabulate(which.max(drop), length(drop)))
    }
    if (debug) {
      dropList <- c(dropList, TipLabels(reduced[[1]])[drop])
      message("Dropping: ", TipLabels(reduced[[1]])[drop],
              " (", which(drop), ")")
    }
    reduced <- keep_and_reduce(tr1, tr2, !drop)
    if (length(reduced) == 1L) {
      reduced <- NULL
    }
    if (debug) {
      if (is.null(reduced[[1]])) {
        plot.new(); plot.new()
      } else {
        plot(reduced[[1]])
        plot(reduced[[2]])
      }
    }
      
    moves <- moves + 1
  }
  
  # Return:
  if (debug) structure(moves, dropList = dropList) else moves
}

# Similar results to phangorn::SPR.dist -- but problem when cutting tree
#' @importFrom TreeTools edge_to_splits
.SPRPairDeOCutter <- function(tree1, tree2, check = TRUE) {
  debug <- isTRUE(getOption("debugSPR", FALSE))
  moves <- 0
  reduced <- ReduceTrees(tree1, tree2, check = check)
  if (debug) {
    dropList <- character(0)
    par(mfrow = 1:2, mai = rep(0.1, 4))
    oldBG <- par(bg = "#eeddcc")
    plot(reduced[[1]])
    nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
    plot(reduced[[2]])
    nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
    par(oldBG)
  }
  
  
  while (!is.null(reduced)) {
    tr1 <- reduced[[1]]
    tr2 <- reduced[[2]]
    edge1 <- tr1[["edge"]]
    edge2 <- tr2[["edge"]]
    labels <- tr1[["tip.label"]]
    nTip <- length(labels)
    sp1 <- edge_to_splits(edge1, PostorderOrder(edge1), labels, nTip = nTip)
    sp2 <- edge_to_splits(edge2, PostorderOrder(edge2), labels, nTip = nTip)

    nSplits <- length(sp1)
    stopifnot(nSplits > 0)
    mmSize <- mismatch_size(sp1, sp2)
    minMismatch <- which.min(mmSize)
    if (mmSize[minMismatch] == 0) {
      agreement <- as.logical(sp1[[.Which1(minMismatch, nSplits)]])
      subtips1 <- agreement
      subtips1[!subtips1][1] <- TRUE
      subtips2 <- !agreement
      subtips2[agreement][1] <- TRUE
      if (debug) {
        action <- paste("Two identical subtrees",
                        names(sp1[[.Which1(minMismatch, nSplits)]]),
                        " = ",
                        names(sp2[[.Which2(minMismatch, nSplits)]])
                        )
        nodelabels("|______", frame = "n", col = "darkred", font = 3,
                   as.integer(names(sp2[[.Which2(minMismatch, nSplits)]])))
        legend("topleft", action, bty = "n")
        message(action)
      }
      
      # The problem with this approach:
      # If subtree one does a clever two-leaf dropper that drops the main tree,
      # we need to subtract one from the overall score.
      
      if (debug) {
        message("> First subtree:")
      }
      submoves1 <- .SPRPairDeOCutter(KeepTipPostorder(reduced[[1]], subtips1),
                            KeepTipPostorder(reduced[[2]], subtips1))
      if (debug) {
        message("> Second subtree:")
      }
      submoves2 <- .SPRPairDeOCutter(KeepTipPostorder(reduced[[1]], subtips2),
                            KeepTipPostorder(reduced[[2]], subtips2))
      return(moves + submoves1 + submoves2)
    }
    split1 <- structure(sp1[.Which1(minMismatch, nSplits), , drop = FALSE],
                        nTip = nTip, class = "Splits")
    split2 <- structure(sp2[.Which2(minMismatch, nSplits), , drop = FALSE],
                        nTip = nTip)
    disagreementSplit <- xor(split1, split2)
    keep <- as.logical(disagreementSplit)
    nKeep <- sum(keep)
    if (nKeep < length(keep) / 2) {
      keep <- !keep
      nKeep <- length(keep) - nKeep
    }

    if (debug) {
      drop <- !keep
      dropList <- c(dropList, labels[drop])
      action <- paste0("Dropping: ", paste(labels[drop], collapse = ", "),
      " (", paste(which(drop), collapse = ", "), ")")
      legend("topleft", action, bty = "n")
      ape::tiplabels("|________", which(drop),
                     frame = "n", col = "darkred", font = 2)
      message(action)
    }
    reduced <- if (nKeep < 4L) {
      NULL
    } else {
      keep_and_reduce(tr1, tr2, keep)
    }
    
    if (length(reduced) == 1L) {
      reduced <- NULL
    }
    
    if (debug) {
      if (is.null(reduced[[1]])) {
        plot.new(); plot.new()
      } else {
        plot(reduced[[1]])
        nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
        plot(reduced[[2]])
        nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
      }
    }
    
    moves <- moves + 1
  }
  
  if (debug) {
    message("> Done.")
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
#' .SPRPairDeO(tree1, tree2)
#' @importFrom TreeTools DropTip TipsInSplits KeepTipPostorder
#' @importFrom TreeTools edge_to_splits
.SPRPairDeO <- function(tree1, tree2, check = TRUE) {
  debug <- isTRUE(getOption("debugSPR", FALSE))
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
      unmatchedSplits <- is.na(matched[["matching"]])
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
