#' Approximate Subtree Prune and Regraft distance
#' 
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
#' 
#' @return `SPRDist()` returns a vector or distance matrix of distances 
#' between trees.
#' 
#' @references \insertAllCited{}
#' 
#' @examples
#' library('TreeTools', quietly = TRUE, warn.conflicts = FALSE)
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
#' @seealso More sophisticated calculation with [\pkg{TBRDist}](
#' https://ms609.github.io/TBRDist/reference/TreeRearrangementDistances.html)
#' functions `USPRDist()` and `ReplugDist()`.
#' 
#' \pkg{phangorn} function \code{\link[phangorn:treedist]{SPR.dist()}} employs
#' the same algorithm but can crash when sent trees of certain formats,
#' and tends to have a longer running time.
#' 
#' @family tree distances
#' @importFrom TreeTools PairwiseDistances Postorder
#' @export
SPRDist <- function (tree1, tree2 = NULL, symmetric) {
  UseMethod("SPRDist")
}

#' @rdname SPRDist
#' @export
SPRDist.phylo <- function (tree1, tree2 = NULL, symmetric) {
  if (is.null(tree2)) {
    NULL
  } else if (inherits(tree2, "phylo")) {
    .SPRPair(tree1, tree2)
  } else {
    vapply(tree2, .SPRPair, double(1), tree1)
  }
}

#' @rdname SPRDist
#' @export
SPRDist.list <- function (tree1, tree2 = NULL, symmetric) {
  if (is.null(tree2)) {
    PairwiseDistances(tree1, .SPRPair)
  } else if (inherits(tree2, 'phylo')) {
    vapply(tree1, .SPRPair, double(1), tree2)
  } else {
    vapply(tree2, SPRDist, double(length(tree1)), tree1)
  }
}

#' @rdname SPRDist
#' @export
SPRDist.multiPhylo <- SPRDist.list

#' @importFrom phangorn SPR.dist
.phangornSPRDist <- function(tree1, tree2 = NULL, symmetric) {
  if (inherits(tree1, 'phylo')) {
    tree1 <- Postorder(tree1)
  } else {
    if (inherits(tree2, 'multiPhylo')) {
      return(vapply(tree2, SPRDist, double(length(tree1)), tree1))
    }
    tree1 <- structure(lapply(tree1, Postorder), class = 'multiPhylo')
  }
  
  if (inherits(tree2, 'phylo')) {
    tree2 <- Postorder(tree2)
  } else if (!is.null(tree2)) {
    tree2 <- structure(lapply(tree2, Postorder), class = 'multiPhylo')
  }
  
  SPR.dist(tree1, tree2)
}



#' @rdname SPRDist
# Using the algorithm of \insertCite{deOliveira2008;textual}{TreeDist}
#' @examples 
#' # de Oliveira Martins et al 2008, fig. 7
#' tree1 <- ape::read.tree(text = "((1, 2), ((a, b), (c, d)), (3, (4, (5, (6, 7)))));")
#' tree2 <- ape::read.tree(text = "((1, 2), 3, (4, (5, (((a, b), (c, d)), (6, 7)))));")
#' plot(tree1)
#' plot(tree2)
#' .SPRPair(tree1, tree2)
#' @importFrom TreeTools DropTip TipsInSplits root_on_node KeepTipFast
#' @export
.SPRConfl <- function(tree1, tree2, debug = FALSE) {
  moves <- 0
  if (debug) dropList <- character(0)
  
  simplified <- Reduce(tree1, tree2)
  if (debug) {
    par(mfrow = 1:2, mai = rep(0.1, 4))
    plot(simplified[[1]])
    plot(simplified[[2]])
    
  }
  
  while (!is.null(simplified)) {
    sp <- as.Splits(simplified)
    nSplits <- length(sp[[1]])
    nTip <- NTip(sp[[1]])
    i <- rep(seq_len(nSplits), nSplits)
    j <- rep(seq_len(nSplits), each = nSplits)
    conf <- confusion(sp[[1]], sp[[2]])
    concave <- colSums(conf == 0)
    
    matches <- concave == 2
    if (any(matches)) {
      agreement <- as.logical(sp[[1]][[i[which.max(matches)]]])
      subtips1 <- agreement
      subtips1[!subtips1][1] <- TRUE
      subtips2 <- !agreement
      subtips2[agreement][1] <- TRUE
      return(moves +
               .SPRPair(KeepTipFast(simplified[[1]], subtips1),
                        KeepTipFast(simplified[[2]], subtips1)) +
               .SPRPair(KeepTipFast(simplified[[1]], subtips2),
                        KeepTipFast(simplified[[2]], subtips2))
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
    
    splitA <- sp[[1]][[i[candidate]]]
    splitB <- sp[[2]][[j[candidate]]]
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
    simplified <- DropTip(simplified, drop)
    simplified <- Reduce(
      root_on_node(simplified[[1]], 1),
      root_on_node(simplified[[2]], 1),
      check = FALSE
    )

    moves <- moves + 1
  }
  
  # Return:
  
  if (debug) list(moves, dropList) else moves
}

# For comparison: not as optimized as phangorn::SPR.dist
#' @importFrom TreeTools xor.Splits
.SPRPairDeO <- function(tree1, tree2) {
  moves <- 0
  simplified <- Reduce(tree1, tree2)
  
  while (!is.null(simplified)) {
    sp <- as.Splits(simplified)
    nSplits <- length(sp[[1]])
    i <- rep(seq_len(nSplits), nSplits)
    j <- rep(seq_len(nSplits), each = nSplits)
    mmSize <- mismatch_size(sp[[1]], sp[[2]])
    minMismatch <- which.min(mmSize)
    if (mmSize[minMismatch] == 0) {
      agreement <- as.logical(sp[[1]][[i[minMismatch]]])
      subtips1 <- agreement
      subtips1[!subtips1][1] <- TRUE
      subtips2 <- !agreement
      subtips2[agreement][1] <- TRUE
      return(moves +
               .SPRPair(KeepTipFast(simplified[[1]], subtips1),
                        KeepTipFast(simplified[[2]], subtips1)) +
               .SPRPair(KeepTipFast(simplified[[1]], subtips2),
                        KeepTipFast(simplified[[2]], subtips2))
             )
    }
    
    disagreementSplit <- xor.Splits(sp[[1]][[i[minMismatch]]],
                                    sp[[2]][[j[minMismatch]]])
    drop <- if (TipsInSplits(disagreementSplit, keep.names = FALSE,
                             smallest = FALSE) != min(mmSize)) {
      !as.logical(disagreementSplit)
    } else {
      as.logical(disagreementSplit)
    }
    # if (sum(drop) > 1) {
    #   drop <- as.logical(tabulate(which.max(drop), length(drop)))
    # }
    simplified <- DropTip(simplified, drop)
    simplified <- Reduce(
      root_on_node(simplified[[1]], 1),
      root_on_node(simplified[[2]], 1),
      check = FALSE
    )
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
  moves
}

.SPRPair <- .SPRConfl
.SPRPair <- .SPRPairDeO
