#' Wrapper for tree distance calculations
#' 
#' Calls tree distance functions from trees or lists of trees
#' 
#' @inheritParams TreeDistance
#' @param Func Tree distance function.
#' @param \dots Additional arguments to `Func`.
#' 
#' @template MRS
#' @keywords internal
#' @importFrom TreeTools as.Splits TipLabels
#' @importFrom utils combn
#' @export
CalculateTreeDistance <- function(Func, tree1, tree2 = NULL,
                                  reportMatching = FALSE, ...) {
  supportedClasses <- c("phylo", "Splits")
  
  single1 <- inherits(tree1, supportedClasses)
  if (!single1 && !inherits(tree1[[1]], supportedClasses)) {
    stop("`tree1` must be a tree, or list of trees, with class ", 
         paste(supportedClasses, collapse = " / "))
  }
  labels1 <- TipLabels(tree1)
  if (is.list(labels1)) {
    if (!all(vapply(labels1[-1], setequal, logical(1), labels1[[1]]))) {
      nTip <- NA
    } else {
      labels1 <- labels1[[1]]
      nTip <- length(labels1)
    }
  } else {
    nTip <- length(labels1)
  }
  
  null2 <- is.null(tree2)
  if (!null2) {
    single2 <- inherits(tree2, supportedClasses)
    if (!single2 && !inherits(tree2[[1]], supportedClasses)) {
      stop("`tree2` must be a tree, or list of trees, with class ", 
           paste(supportedClasses, collapse = " / "))
    }
    labels2 <- TipLabels(tree2)
    if (is.list(labels2)) {
      if (!all(vapply(labels2[-1], setequal, logical(1), labels2[[1]]))) {
        nTip <- NA
      } else {
        labels2 <- labels2[[1]]
      }
    }
    
    if (!setequal(labels1, labels2)) {
      nTip <- NA
    }
  }
  
  if (single1) {
    if (null2) {
      dist(0)
    } else if (single2) {
      .SplitDistanceOneOne(Func, tree1, tree2, tipLabels = labels1, nTip,
                           reportMatching = reportMatching, ...)
    } else {
      .SplitDistanceOneMany(Func, oneSplit = tree1, manySplits = tree2,
                            tipLabels = labels1, nTip = nTip, ...)
    }
  } else {
    if (is.null(tree2)) {
      .SplitDistanceAllPairs(Func, tree1, tipLabels = labels1, nTip, ...)
    } else if (single2) {
      .SplitDistanceOneMany(Func, oneSplit = tree2, manySplits = tree1,
                            tipLabels = labels2, nTip = nTip, ...)
    } else {
      .SplitDistanceManyMany(Func, tree1, tree2, tipLabels = labels1,
                             nTip = nTip, ...)
    }
  }
}

.SplitDistanceOneOne <- function(Func, split1, split2, tipLabels, 
                                 nTip = length(tipLabels), reportMatching,
                                 ...) {
  if (is.na(nTip)) {
    common <- intersect(TipLabels(split1), TipLabels(split2))
    split1 <- KeepTip(split1, common)
    split2 <- KeepTip(split2, common)
    nTip <- length(common)
  }
  if (nTip < 4) {
    0
  } else {
    Func(as.Splits(split1, asSplits = reportMatching),
         as.Splits(split2, tipLabels = tipLabels, asSplits = reportMatching),
         nTip = nTip, reportMatching = reportMatching, ...)
  }
}

#' @importFrom stats setNames
.SplitDistanceOneMany <- function(Func, oneSplit, manySplits, 
                                  tipLabels, nTip = length(tipLabels), ...) {
  if (is.na(nTip)) {
    oneLabels <- TipLabels(oneSplit)
    manyLabels <- TipLabels(manySplits)
    common <- if (is.list(manyLabels)) {
      lapply(manyLabels, intersect, oneLabels)
    } else {
      list(intersect(manyLabels, oneLabels))
    }
    setNames(unlist(.mapply(function(s2, keep) {
      Func(
        as.Splits(KeepTip(oneSplit, keep), tipLabels = keep, asSplits = FALSE),
        as.Splits(KeepTip(s2, keep), tipLabels = keep, asSplits = FALSE),
        nTip = length(keep),
        ...
      )
    }, dots = list(manySplits, common), MoreArgs = NULL)), names(manySplits))
  } else {
    s1 <- as.Splits(oneSplit, tipLabels = tipLabels, asSplits = FALSE)
    vapply(manySplits,
           function(s2) Func(s1, as.Splits(s2, tipLabels = tipLabels,
                                           asSplits = FALSE),
                              nTip = nTip, ...),
           double(1))
  }
}

#' @importFrom parallel parCapply
#' @importFrom cli cli_progress_bar cli_progress_update
.SplitDistanceAllPairs <- function(Func, splits1, tipLabels,
                                   nTip = length(tipLabels), ...) {
  cluster <- getOption("TreeDist-cluster")

  if (is.na(nTip)) {
    splits <- lapply(splits1, as.Splits)
    
    .PairDist <- function(i) {
      s1 <- splits[[i[[1]]]]
      s2 <- splits[[i[[2]]]]
      common <- intersect(TipLabels(s1), TipLabels(s2))
      Func(KeepTip(s1, common), KeepTip(s2, common),
           nTip = length(common), reportMatching = FALSE, ...)
    }
    
  } else {
    splits <- as.Splits(splits1, tipLabels = tipLabels, asSplits = FALSE)
    
    .PairDist <- function(i) {
      Func(splits[[i[[1]]]], splits[[i[[2]]]],
           nTip = nTip, reportMatching = FALSE, ...)
    }
    
  }
  
  nSplits <- length(splits)
  is <- combn(seq_len(nSplits), 2)
  
  if (is.null(cluster)) {
    cli_progress_bar("Calculating distances", total = ncol(is))
  }
  
  .CliPairDist <- function(i) { # TODO if cli possible from parallel, merge.
    cli_progress_update(1, .envir = parent.frame(2))
    .PairDist(i)
  }
  
  ret <- structure(class = "dist", Size = nSplits,
                   Labels = names(splits1),
                   Diag = FALSE, Upper = FALSE,
                   if (is.null(cluster)) {
                     apply(is, 2, .CliPairDist)
                   } else {
                     parCapply(cluster, is, .PairDist)
                   })
  # Return:
  ret
}

#' @importFrom stats setNames
.SplitDistanceManyMany <- function(Func, splits1, splits2, 
                                   tipLabels, nTip = length(tipLabels), ...) {
  
  if (is.na(nTip)) {
    tipLabels <- union(unlist(tipLabels, use.names = FALSE),
                       unlist(TipLabels(splits2), use.names = FALSE))
    splits1 <- as.Splits(splits1, tipLabels = tipLabels, asSplits = TRUE)
    splits2 <- as.Splits(splits2, tipLabels = tipLabels, asSplits = TRUE)
    vapply(splits1, function(s1) {
        l1 <- TipLabels(s1)
        vapply(splits2, function(s2) {
          common <- intersect(l1, TipLabels(s2))
          Func(KeepTip(s1, common), KeepTip(s2, common),
               nTip = length(common))#, ...)
        }, double(1))
      },
      setNames(double(length(splits2)), names(splits2))
    )
  } else {
    splits1 <- as.Splits(splits1, tipLabels = tipLabels, asSplits = FALSE)
    splits2 <- as.Splits(splits2, tipLabels = tipLabels, asSplits = FALSE)
    matrix(
      unlist(.mapply(Func, list(rep(splits2, each = length(splits1)),
                         splits1, nTip = nTip, ...), NULL)),
      length(splits1),
      length(splits2),
      dimnames = list(names(splits1), names(splits2))
    )
  }
}

#' Calculate distance between trees, or lists of trees
#' @template MRS
#' @importFrom TreeTools TipLabels
#' @param checks Logical specifying whether to perform basic sanity checks to
#' avoid crashes in C++.
#' @keywords internal
#' @seealso [`CalculateTreeDistance`]
#' @export
.TreeDistance <- function(Func, tree1, tree2, checks = TRUE, ...) {
  single1 <- inherits(tree1, "phylo")
  labels1 <- TipLabels(tree1, single = TRUE)
  nTip <- length(labels1)
  
  single2 <- inherits(tree2, "phylo")
  labels2 <- TipLabels(tree2, single = TRUE)
  
  if (checks) {
    if (length(setdiff(labels1, labels2)) > 0) {
      stop("Leaves must bear identical labels.")
    }
    
    if (!single1) {
      .CheckLabelsSame(lapply(tree1, TipLabels))
    }
    
    if (!single2) {
      .CheckLabelsSame(lapply(tree2, TipLabels))
    }
  }
  
  if (single1) {
    if (single2) {
      Func(tree1, tree2, tipLabels = labels1, nTip = nTip, ...)
    } else {
      .TreeDistanceOneMany(Func, oneTree = tree1, manyTrees = tree2,
                           tipLabels = labels1, nTip = nTip, ...)
    }
  } else {
    if (single2) {
      .TreeDistanceOneMany(Func, oneTree = tree2, tipLabels = labels2, 
                           nTip = nTip, manyTrees = tree1, ...)
    } else {
      .TreeDistanceManyMany(Func, tree1, tree2, tipLabels = labels1,
                            nTip = nTip, ...)
    }
  }
}

.TreeDistanceOneMany <- function(Func, oneTree, manyTrees, tipLabels, 
                                  nTip = length(tipLabels), 
                                  FUN.VALUE = Func(manyTrees[[1]], oneTree,
                                                   nTip = nTip, 
                                                   tipLabels = tipLabels, ...),
                                  ...) {
  vapply(manyTrees, Func, oneTree, tipLabels = tipLabels,
                            nTip = nTip, ..., FUN.VALUE = FUN.VALUE)
}

.TreeDistanceManyMany <- function(Func, trees1, trees2, tipLabels, 
                                   nTip = length(tipLabels),
                                   FUN.VALUE = Func(trees1[[1]], trees2[[1]],
                                                    nTip = nTip, 
                                                    tipLabels = tipLabels, ...),
                                   ...) {
  if (identical(trees1, trees2)) {
    CompareAll(trees1, Func, nTip = nTip, tipLabels = tipLabels,
               FUN.VALUE = FUN.VALUE, ...)
  } else {
    seqAlong1 <- seq_along(trees1)
    names(seqAlong1) <- names(trees1)
    value <- vapply(seqAlong1, function(x) FUN.VALUE, FUN.VALUE)
    ret <- vapply(trees2, .TreeDistanceOneMany, Func = Func, manyTrees = trees1,
           tipLabels = tipLabels, nTip = nTip, FUN.VALUE = value, ...)
    if (length(dim(ret)) == 3) {
      aperm(ret, c(2, 3, 1))
    } else {
      ret
    }
  }
}

.CheckLabelsSame <- function(labelList) {
  nTip <- unique(lengths(labelList))
  if (length(nTip) != 1) {
    stop("All trees must contain the same number of leaves.")
  }
  tipLabel <- unique(lapply(labelList, sort))
  if (length(tipLabel) != 1L) {
    stop("All trees must bear identical labels. Found:\n >  ", 
         paste0(lapply(tipLabel, paste, collapse = " "), 
               " (", lapply(tipLabel, class), ")",
                collapse = "\n >  "))
  }
}

#' Entropy in bits
#' 
#' Calculate the entropy of a vector of probabilities, in bits.
#' Probabilities should sum to one.
#' Probabilities equalling zero will be ignored.
#' 
#' @param \dots Numerics or numeric vector specifying probabilities of outcomes.
#' 
#' @return `Entropy()` returns the entropy of the specified probabilities, 
#' in bits.
#' 
#' @examples
#' Entropy(1/2, 0, 1/2) # = 1
#' Entropy(rep(1/4, 4)) # = 2
#' @template MRS
#' @export
Entropy <- function(...) {
  p <- c(...)
  p <- p[p > 0]
  -sum(p * log2(p))
}

#' Distances between each pair of trees
#' 
#' Calculate the distance between each tree in a list, and each other tree
#' in the same list.
#' 
#' `CompareAll()` is not limited to tree comparisons:
#' `Func` can be any symmetric function.
#'
#' @param x List of trees, in the format expected by `Func()`.
#' @param Func distance function returning distance between two trees,
#' e.g. [`path.dist()`][phangorn::treedist].
#' @param FUN.VALUE Format of output of `Func()`, to be passed to [`vapply()`]. 
#' If unspecified, calculated by running `Func(x[[1]], x[[1]])`.
#' @param \dots Additional parameters to pass to `Func()`.
#' @return `CompareAll()` returns a distance matrix of class `dist` detailing
#' the distance between each pair of trees.
#' Identical trees are assumed to have zero distance.
#' 
#' @examples
#' # Generate a list of trees to compare
#' library("TreeTools", quietly = TRUE)
#' trees <- list(bal1 = BalancedTree(1:8), 
#'               pec1 = PectinateTree(1:8),
#'               pec2 = PectinateTree(c(4:1, 5:8)))
#' 
#' # Compare each tree with each other tree
#' CompareAll(trees, NNIDist)
#'   
#' # Providing FUN.VALUE yields a modest speed gain:
#' dist <- CompareAll(trees, NNIDist, FUN.VALUE = integer(7))
#'   
#' # View distances as a matrix
#' as.matrix(dist$lower)
#' @template MRS
#' @family pairwise tree distances
#' @importFrom cli cli_progress_bar cli_progress_update
#' @importFrom parallel parLapply
#' @importFrom stats dist
#' @export
CompareAll <- function(x, Func, FUN.VALUE = Func(x[[1]], x[[1]], ...),
                        ...) {
  nTree <- length(x)
  countUp <- seq_len(nTree - 1)
  i <- x[rep(countUp, rev(countUp))]
  j <- x[unlist(sapply(countUp, function(n) n + seq_len(nTree - n)))]
  
  cluster <- getOption("TreeDist-cluster")
  ret <- if (is.null(cluster)) {
    cli_progress_bar("Comparing", total = length(i))
    vapply(seq_along(i), function(k) {
      cli_progress_update(1, .envir = parent.frame(2))
      Func(i[[k]], j[[k]], ...)
      }, FUN.VALUE)
  } else {
    do.call("cbind", 
            parLapply(cluster, seq_along(i), 
                      function(k, Func) Func(i[[k]], j[[k]]),
                      Func = Func))
  }
  
  .WrapReturn <- function(dists) {
    structure(dists,
              Size = nTree,
              Labels = names(x),
              Diag = FALSE,
              Upper = TRUE,
              class = "dist")
  }
  
  # Return:
  if (length(FUN.VALUE) == 1) {
    .WrapReturn(ret)
  } else {
    structure(lapply(seq_len(dim(ret)[[1]]), 
                     function(i) .WrapReturn(unlist(ret[i, ]))),
              names = rownames(ret))
  }
}

#' Normalize tree distances
#' 
#' `NormalizeInfo()` is an internal function used to normalize information
#' against a reference, such as the total information present in a pair of
#' trees.
#' 
#' The unnormalized value(s) are normalized by dividing by a denominator
#' calculated based on the `how` parameter.  Valid options include:
#' 
#' \describe{
#'  \item{`FALSE`}{No normalization is performed; the unnormalized values
#'  are returned.}
#'  \item{`TRUE`}{Unless `infoInBoth` is specified, the information in
#'  each tree is computed using `InfoInTree()`, and the two values combined
#'  using `Combine()`.}
#'  \item{A numeric value, vector or matrix}{`how` is used as the denominator;
#'  the returned value is `unnormalized / how`.}
#'  \item{A function}{Unless `infoInBoth` is specified, the information in
#'  each tree is computed using `InfoInTree()`, and the two values combined
#'  using `how`.  `NormalizeInfo(how = Func)` is thus equivalent to
#'  `NormalizeInfo(how = TRUE, Combine = Func)`.}
#' }
#' 
#' @param unnormalized Numeric value, vector or matrix to be normalized.
#' @param tree1,tree2 Trees from which `unnormalized` was calculated.
#' @param InfoInTree Function to calculate the information content of each tree.
#' @param infoInBoth Optional numeric specifying information content of both
#' trees independently. If unspecified (`NULL`), this will be calculated using
#' the method specified by `how`.
#' @param how Method for normalization, perhaps specified using the `normalize`
#' argument to a tree distance function.  See details for options.
#' @param \dots Additional parameters to `InfoInTree()` or `how()`.
#' @returns `NormalizeInfo()` returns an object corresponding to the normalized
#' values of `unnormalized`.
#' @examples
#' library("TreeTools", quietly = TRUE)
#' pair1 <- c(BalancedTree(9), StarTree(9))
#' pair2 <- c(BalancedTree(9), PectinateTree(9))
#' 
#' # We'll let the number of nodes define the total information in a tree
#' Nnode(pair1)
#' Nnode(pair2)
#' 
#' # Let's normalize a unit distance
#' rawDist <- cbind(c(1, 1), c(1, 1))
#' 
#' # With `Combine = "+"`, the maximum distance is the sum of
#' # the information in each tree
#' denominator <- outer(Nnode(pair1), Nnode(pair2), "+")
#' 
#' NormalizeInfo(rawDist, pair1, pair2, InfoInTree = ape::Nnode, Combine = "+")
#' rawDist / denominator
#' 
#' 
#' # A denominator can be specified manually using `how`:
#' NormalizeInfo(rawDist, pair1, pair2, InfoInTree = ape::Nnode, how = 16)
#' rawDist / 16
#' 
#' 
#' # `how` also allows the denominator to be computed from trees:
#' outer(Nnode(pair1), Nnode(pair2), pmin)
#' NormalizeInfo(rawDist, pair1, pair2, InfoInTree = ape::Nnode, how = pmin)
#' rawDist / outer(Nnode(pair1), Nnode(pair2), pmin)
#' 
#' @keywords internal
#' @template MRS
#' @importFrom TreeTools KeepTip TipLabels
#' @export
NormalizeInfo <- function(unnormalized, tree1, tree2, InfoInTree,
                          infoInBoth = NULL, how = TRUE, Combine = "+", ...) {
  
  CombineInfo <- function(tree1Info, tree2Info, Combiner = Combine,
                          pairwise = FALSE) {
    if (length(tree1Info) == 1 || length(tree2Info) == 1 || pairwise) {
      # TODO When requriring R4.0, remove match.fun - which is now part of
      # .mapply
      unlist(.mapply(match.fun(Combiner),
                     dots = list(tree1Info, tree2Info), NULL))
    } else {
      ret <- outer(tree1Info, tree2Info, Combiner)
      if (inherits(unnormalized, "dist")) ret[lower.tri(ret)] else ret
    }
  }
  
  lab1 <- TipLabels(tree1)
  lab2 <- TipLabels(tree2)
  sameLabels <- .AllTipsSame(lab1, lab2)
  
  if (!sameLabels) {
    trees <- .SharedOnly(tree1, tree2, lab1, lab2)
    tree1 <- trees[[1]]
    tree2 <- trees[[2]]
  }
  
  if (is.logical(how)) {
    if (how == FALSE) {
      return(unnormalized)
    } else {
      if (is.null(infoInBoth)) {
        info1 <- InfoInTree(tree1, ...)
        info2 <- if (is.null(tree2)) {
          info1
        } else {
          InfoInTree(tree2, ...)
        }
        infoInBoth <- CombineInfo(info1, info2, pairwise = !sameLabels)
      }
    }
  } else if (is.function(how)) {
    if (is.null(infoInBoth)) {
      info1 <- InfoInTree(tree1, ...)
      info2 <- if (is.null(tree2)) info1 else InfoInTree(tree2, ...)
      infoInBoth <- CombineInfo(info1, info2, Combiner = how,
                                pairwise = !sameLabels)
    }
  } else {
    infoInBoth <- how
  }
  # Return:
  unnormalized / infoInBoth
}

# We only call this function when not all trees contain identical leaf sets
#' @importFrom TreeTools KeepTip TipLabels
.SharedOnly <- function(tree1, tree2,
                        lab1 = TipLabels(tree1),
                        lab2 = TipLabels(tree2)) {
  if (is.null(tree2)) {
    # Case: N trees vs themselves
    # We require a triangular matrix suitable for a dist object
    if (.MultipleTrees(tree1)) {
      pairs <- combn(seq_along(tree1), 2)
      nPairs <- dim(pairs)[[2]]
      ret <- list(vector("list", nPairs), vector("list", nPairs))
      
      for (n in seq_len(nPairs)) {
        i <- pairs[1, n]
        j <- pairs[2, n]
        common <- intersect(lab1[[i]], lab1[[j]])
        ret[[1]][[n]] <- KeepTip(tree1[[i]], common)
        ret[[2]][[n]] <- KeepTip(tree1[[j]], common)
      }
      
      # Return:
      ret
    } else {
      return(list(NULL, NULL))
    }
  } else if (.MultipleTrees(tree1)) {
    if (.MultipleTrees(tree2)) {
      # Case: N vs M
      n1 <- length(tree1)
      n2 <- length(tree2)
      nPairs <- n1 * n2
      ret <- list(vector("list", nPairs), vector("list", nPairs))
      pair1 <- rep(tree1, each = n2)
      pair2 <- rep.int(tree2, times = n1)
      lab1 <- if (is.list(lab1)) {
        rep(lab1, each = n2)
      } else {
        rep.int(list(lab1), times = n1 * n2)
      }
      lab2 <- if (is.list(lab2)) {
        rep.int(lab2, times = n1)
      } else {
        rep.int(list(lab2), times = n1 * n2)
      }
      for (i in seq_len(nPairs)) {
        common <- intersect(lab1[[i]], lab2[[i]])
        ret[[1]][[i]] <- KeepTip(pair1[[i]], common)
        ret[[2]][[i]] <- KeepTip(pair2[[i]], common)
      }
      # Return:
      ret
    } else {
      # Case: N vs 1
      if (inherits(tree2, "multiPhylo")) {
        tree2 <- tree2[[1]]
      }
      commonLeaves <- lapply(if (is.list(lab1)) lab1 else list(lab1),
                             intersect, lab2)
      list(.mapply(KeepTip, list(tree1, commonLeaves), NULL),
           lapply(commonLeaves, function (common) KeepTip(tree2, common)))
    }
  } else {
    if (.MultipleTrees(tree2)) {
      # Case: 1 vs N
      if (inherits(tree1, "multiPhylo")) {
        tree1 <- tree1[[1]]
      }
      commonLeaves <- lapply(if (is.list(lab2)) lab2 else list(lab2),
                             intersect, lab1)
      list(lapply(commonLeaves, function (common) KeepTip(tree1, common)),
           .mapply(KeepTip, list(tree2, commonLeaves), NULL))
    } else {
      # Case: 1 vs 1
      commonLeaves <- intersect(lab1, lab2)
      list(KeepTip(tree1, commonLeaves), KeepTip(tree2, commonLeaves))
    }
  }
}

.MultipleTrees <- function(tree) {
  is.list(tree) && inherits(tree[[1]], "phylo") && length(tree) > 1
}

#' List clades as text
#' @param splits1,splits2 Logical matrices with columns specifying membership
#' of each corresponding matched clade.
#' @return `ReportMatching` returns a character vector describing each pairing 
#' in a matching.
#'   
#' @seealso [`VisualizeMatching`]
#' @template MRS
#' @keywords internal
#' @export
ReportMatching <- function(splits1, splits2, realMatch = TRUE) {
  paste(as.character(splits1), ifelse(realMatch, "=>", ".."),
        as.character(splits2))
}
