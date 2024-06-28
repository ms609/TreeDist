#' Approximate the Subtree Prune and Regraft (SPR) distance.
#' 
#' `SPRDist()` calculates an upper bound on the SPR distance between trees
#' using the heuristic method of \insertCite{deOliveira2008;textual}{TreeDist}.
#' Other approximations are available
#' \insertCite{@e.g. @Hickey2008, @Goloboff2008SPR, @Whidden2018}{TreeDist}.
#' 
#' @template tree12ListParams
#' @param symmetric Ignored (redundant after fix of
#' [phangorn#97](https://github.com/KlausVigo/phangorn/issues/97)).
#' @param method Character specifying which method to use to approximate the
#' SPR distance.  Currently defaults to "deOliveira", the only accepted option;
#' a new method will be available soon.
#' 
#' @return `SPRDist()` returns a vector or distance matrix of distances 
#' between trees.
#' 
#' @references \insertAllCited{}
#' 
#' @examples
#' library("TreeTools", quietly = TRUE)
#' 
#' SPRDist(BalancedTree(7), PectinateTree(7))
#' 
#' SPRDist(BalancedTree(7), as.phylo(0:2, 7))
#' SPRDist(as.phylo(0:2, 7), PectinateTree(7))
#'
#' SPRDist(list(bal = BalancedTree(7), pec = PectinateTree(7)),
#'         as.phylo(0:2, 7))
#'
#' CompareAll(as.phylo(30:33, 8), SPRDist)
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
SPRDist <- function (tree1, tree2 = NULL, symmetric, method = "deOliveira") {
  UseMethod("SPRDist")
}

#' @rdname SPRDist
#' @export
SPRDist.phylo <- function (tree1, tree2 = NULL, symmetric, method = "deOliveira") {
  if (is.null(tree2)) {
    NULL
  } else if (inherits(tree2, "phylo")) {
    .SPRFunc(method)(tree1, tree2)
  } else {
    vapply(tree2, .SPRFunc(method), double(1), tree1)
  }
}

.SPRFunc <- function(method) {
  switch(pmatch(tolower(method), c("deoliveira", "confl", "experiment")),
         .SPRPairDeOCutter, .SPRConfl, .SPRExperiment)
}

#' @rdname SPRDist
#' @export
SPRDist.list <- function (tree1, tree2 = NULL, symmetric, method = "deOliveira") {
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


.SPRExperiment <- function(tree1, tree2, check = TRUE, debug = FALSE) {
  moves <- 0
  if (debug) dropList <- character(0)
  
  simplified <- Reduce(tree1, tree2, check = check)
  if (debug) {
    dropList <- character(0)
    par(mfrow = 1:2, mai = rep(0.1, 4))
    oldBG <- par(bg = "#eeddcc")
    plot(simplified[[1]])
    nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
    plot(simplified[[2]])
    nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
    par(oldBG)
  }
  
  while (!is.null(simplified)) {
    tr1 <- simplified[[1]]
    tr2 <- simplified[[2]]
    edge1 <- tr1[["edge"]]
    edge2 <- tr2[["edge"]]
    labels <- tr1[["tip.label"]]
    nTip <- length(labels)
    sp1 <- edge_to_splits(edge1, PostorderOrder(edge1), labels, nTip = nTip)
    sp2 <- edge_to_splits(edge2, PostorderOrder(edge2), labels, nTip = nTip)
    nSplits <- length(sp1)
    
    conf <- confusion(sp1, sp2)
    if (debug) {
      dimnames(conf) <- list(
        c("ab", "aB", "Ab", "AB"),
        names(sp1),
        names(sp2))
    }
    concave <- colSums(conf == 0)
    
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
      return(moves +
               .SPRPair(KeepTipPostorder(tr1, subtips1),
                        KeepTipPostorder(tr2, subtips1)) +
               .SPRPair(KeepTipPostorder(tr1, subtips2),
                        KeepTipPostorder(tr2, subtips2))
      )
    }
    .Is1 <- function (i, j) {
      hitHere <- logical(attr(sp1, "nTip"))
      if (conf[1, i, j] == 1) {
        hitHere <- hitHere | as.logical(sp1[[i]] & sp2[[j]])
      }
      if (conf[2, i, j] == 1) {
        hitHere <- hitHere | as.logical(sp1[[i]] & !sp2[[j]])
      }
      if (conf[3, i, j] == 1) {
        hitHere <- hitHere | as.logical(!sp1[[i]] & sp2[[j]])
      }
      if (conf[4, i, j] == 1) {
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
      if (conf[1, i, j] > 1) {
        as.logical(!sp1[[i]] & !sp2[[j]])
      } else 
      if (conf[2, i, j] > 1) {
        as.logical(!sp1[[i]] & sp2[[j]])
      } else
      if (conf[3, i, j] > 1) {
        as.logical(sp1[[i]] & !sp2[[j]])
      } else
      if (conf[4, i, j] > 1) {
        as.logical(sp1[[i]] & sp2[[j]])
      })
    }
    
    nits <- which(apply(conf, 2:3, function (x) sum(0:2 %in% x)) == 3)
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
      twits <- which(apply(conf, 2:3, function (x) sum(x == 1) > 1))
      if (length(twits)) {
        twitDrops <- unlist(sapply(twits, .FindDrops))
        keep <- !tabulate(which.max(tabulate(twitDrops)), nTip)
      
        # flits <- which(apply(conf, 2:3, function (x) sum(x == 1) == 3))
        # flitDrops <- vapply(flits, .FindOverlap, integer(1))
        # nitDrops <- unique(flitDrops)
        if (debug) {
          message("    Flit candidates: ",
                  paste(labels[nitDrops], collapse = ", "),
                  "; duplicates: ",
                  paste(labels[flitDrops[duplicated(as.integer(flitDrops))]],
                        collapse = ", ")
          )
        }
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
    simplified <- if (nKeep < 4L) {
      NULL
    } else {
      keep_and_reduce(tr1, tr2, keep)
    }
    if (length(simplified) == 1) {
      simplified <- NULL
    }
    
    if (debug) {
      if (is.null(simplified)) {
        plot.new(); plot.new()
      } else {
        plot(simplified[[1]])
        nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
        plot(simplified[[2]])
        nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
      }
    }
    
    moves <- moves + 1
  }
  
  # Return:
  
  if (debug) list(moves, dropList) else moves
}


.SPRConfl <- function(tree1, tree2, check = TRUE, debug = FALSE) {
  moves <- 0
  if (debug) dropList <- character(0)
  
  simplified <- Reduce(tree1, tree2, check = check)
  if (debug) {
    par(mfrow = 1:2, mai = rep(0.1, 4))
    plot(simplified[[1]])
    plot(simplified[[2]])
  }
  
  while (!is.null(simplified)) {
    tr1 <- simplified[[1]]
    tr2 <- simplified[[2]]
    edge1 <- tr1[["edge"]]
    edge2 <- tr2[["edge"]]
    labels <- tr1[["tip.label"]]
    nTip <- length(labels)
    sp1 <- edge_to_splits(edge1, PostorderOrder(edge1), labels, nTip = nTip)
    sp2 <- edge_to_splits(edge2, PostorderOrder(edge2), labels, nTip = nTip)
    nSplits <- length(sp1)
    
    conf <- confusion(sp1, sp2)
    concave <- colSums(conf == 0)
    
    matches <- concave == 2
    if (any(matches)) {
      agreement <- as.logical(sp1[[.Which1(which.max(matches), nSplits)]])
      subtips1 <- agreement
      subtips1[!subtips1][1] <- TRUE
      subtips2 <- !agreement
      subtips2[agreement][1] <- TRUE
      return(moves +
               .SPRPair(KeepTipPostorder(tr1, subtips1),
                        KeepTipPostorder(tr2, subtips1)) +
               .SPRPair(KeepTipPostorder(tr1, subtips2),
                        KeepTipPostorder(tr2, subtips2))
             )
    }
    
    confInf <- conf
    confInf[conf == 0] <- Inf
    confMin <- apply(confInf, 2:3, min)
    minConf <- min(confMin[confMin > 0])
    if (debug && minConf > 1) {
      message("Minimum conflict: ", minConf)
    }
    h <- apply(conf / nTip, 2:3, Entropy)
    minH <- min(h[confMin == minConf])
    maxH <- max(h[confMin == minConf])
    
    candidate <- which.max(h == maxH)
    
    splitA <- sp1[[.Which1(candidate, nSplits)]]
    splitB <- sp2[[.Which2(candidate, nSplits)]]
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
      dropList <- c(dropList, TipLabels(simplified[[1]])[drop])
      message("Dropping: ", TipLabels(simplified[[1]])[drop],
              " (", which(drop), ")")
    }
    simplified <- keep_and_reduce(tr1, tr2, !drop)
    if (length(simplified) == 1L) {
      simplified <- NULL
    }
    if (debug) {
      if (is.null(simplified[[1]])) {
        plot.new(); plot.new()
      } else {
        plot(simplified[[1]])
        plot(simplified[[2]])
      }
    }
      
    moves <- moves + 1
  }
  
  # Return:
  
  if (debug) list(moves, dropList) else moves
}

# Similar results to phangorn::SPR.dist -- but problem when cutting tree
#' @importFrom TreeTools edge_to_splits
.SPRPairDeOCutter <- function(tree1, tree2, check = TRUE, debug = FALSE) {
  moves <- 0
  simplified <- Reduce(tree1, tree2, check = check)
  if (debug) {
    dropList <- character(0)
    par(mfrow = 1:2, mai = rep(0.1, 4))
    oldBG <- par(bg = "#eeddcc")
    plot(simplified[[1]])
    nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
    plot(simplified[[2]])
    nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
    par(oldBG)
  }
  
  
  while (!is.null(simplified)) {
    tr1 <- simplified[[1]]
    tr2 <- simplified[[2]]
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
      submoves1 <- .SPRPairDeOCutter(KeepTipPostorder(simplified[[1]], subtips1),
                            KeepTipPostorder(simplified[[2]], subtips1),
                            debug = debug)
      if (debug) {
        message("> Second subtree:")
      }
      submoves2 <- .SPRPairDeOCutter(KeepTipPostorder(simplified[[1]], subtips2),
                            KeepTipPostorder(simplified[[2]], subtips2),
                            debug = debug)
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
    simplified <- if (nKeep < 4L) {
      NULL
    } else {
      keep_and_reduce(tr1, tr2, keep)
    }
    
    if (length(simplified) == 1L) {
      simplified <- NULL
    }
    
    if (debug) {
      if (is.null(simplified[[1]])) {
        plot.new(); plot.new()
      } else {
        plot(simplified[[1]])
        nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
        plot(simplified[[2]])
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

# An attempt to reproduce the phangorn results using the actual algorithm
# described, in which matched edges are not considered further.
# Using the algorithm of \insertCite{deOliveira2008;textual}{TreeDist}
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
.SPRPairDeO <- function(tree1, tree2, check = TRUE, debug = FALSE) {
  moves <- 0
  simplified <- Reduce(tree1, tree2, check = check)
  if (debug) {
    dropList <- character(0)
    par(mfrow = 1:2, mai = rep(0.1, 4))
    oldBG <- par(bg = "#eeddcc")
    plot(simplified[[1]])
    nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
    plot(simplified[[2]])
    nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
    par(oldBG)
  }
  
  
  while (!is.null(simplified)) {
    tr1 <- simplified[[1]]
    tr2 <- simplified[[2]]
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
    stopifnot(cpp_robinson_foulds_distance(sp1, sp2, nTip)$score ==
                matched$score)
    
    nSplits <- length(sp1)
    stopifnot(nSplits > 0)
    mmSize <- mismatch_size(sp1, sp2)
    sapply(which(mmSize == 0), .Which1, nSplits)
    sapply(which(mmSize == 0), .Which2, nSplits)
    minMismatch <- which.min(mmSize)
    stopifnot(mmSize[minMismatch] > 0)
    
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
                     frame = "n", col = "red", font = 2)
      message(action)
    }
    simplified <- if (nKeep < 4L) {
      NULL
    } else {
      keep_and_reduce(tr1, tr2, keep)
    }
    
    if (length(simplified) == 1L) {
      simplified <- NULL
    }
    
    if (debug) {
      if (is.null(simplified[[1]])) {
        plot.new(); plot.new()
      } else {
        plot(simplified[[1]])
        nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
        plot(simplified[[2]])
        nodelabels(cex = 0.7, bg = NA, frame = "n", adj = 1.5)
      }
    }
    
    moves <- moves + 1
  }
  
  # Return:
  moves
}
