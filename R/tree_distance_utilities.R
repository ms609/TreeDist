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
CalculateTreeDistance <- function (Func, tree1, tree2 = NULL,
                                   reportMatching = FALSE, ...) {
  supportedClasses <- c('phylo', 'Splits')
  
  single1 <- inherits(tree1, supportedClasses)
  if (!single1 && !inherits(tree1[[1]], supportedClasses)) {
    stop("`tree1` must be a tree, or list of trees, with class ", 
         paste(supportedClasses, collapse = ' / '))
  }
  labels1 <- TipLabels(tree1, single = TRUE)
  nTip <- length(labels1)
  
  null2 <- is.null(tree2)
  if (!null2) {
    single2 <- inherits(tree2, supportedClasses)
    if (!single2 && !inherits(tree2[[1]], supportedClasses)) {
      stop("`tree2` must be a tree, or list of trees, with class ", 
           paste(supportedClasses, collapse = ' / '))
    }
    labels2 <- TipLabels(tree2, single = TRUE)
    
    if (length(setdiff(labels1, labels2)) > 0) {
      stop("Leaves must bear identical labels.")
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
                            tipLabels = labels2, ...)
    } else {
      .SplitDistanceManyMany(Func, tree1, tree2,
                             tipLabels = labels1, ...)
    }
  }
}

.SplitDistanceOneOne <- function (Func, split1, split2, tipLabels, 
                                 nTip = length(tipLabels), reportMatching, ...) {
  Func(as.Splits(split1, asSplits = reportMatching),
       as.Splits(split2, tipLabels = tipLabels, asSplits = reportMatching),
       nTip = nTip, reportMatching = reportMatching, ...)
}

.SplitDistanceOneMany <- function (Func, oneSplit, manySplits, 
                                  tipLabels, nTip = length(tipLabels), ...) {
  s1 <- as.Splits(oneSplit, tipLabels = tipLabels, asSplits = FALSE)
  vapply(manySplits,
         function (s2) Func(s1, as.Splits(s2, tipLabels = tipLabels,
                                          asSplits = FALSE),
                            nTip = nTip, ...),
         double(1))
}

.SplitDistanceAllPairs <- function (Func, splits1, tipLabels,
                                    nTip = length(tipLabels), ...) {
  splits <- as.Splits(splits1, tipLabels = tipLabels, asSplits = FALSE)
  nSplits <- length(splits)
  is <- combn(seq_len(nSplits), 2)
  
  ret <- structure(class = 'dist', Size = nSplits,
                   diag = FALSE, upper = FALSE,
                   apply(is, 2, function (i)
                     Func(splits[[i[1]]], splits[[i[2]]],
                          nTip = nTip, reportMatching = FALSE, ...)))
  # Return:
  ret
}

.SplitDistanceManyMany <- function (Func, splits1, splits2, 
                                    tipLabels, nTip = length(tipLabels), ...) {
  splits1 <- as.Splits(splits1, tipLabels = tipLabels, asSplits = FALSE)
  splits2 <- as.Splits(splits2, tipLabels = tipLabels, asSplits = FALSE)
  matrix(mapply(Func, rep(splits2, each = length(splits1)), splits1,
                nTip = nTip, ...),
         length(splits1), length(splits2),
         dimnames = list(names(splits1), names(splits2)))
}

#' Calculate distance between trees, or lists of trees
#' @template MRS
#' @importFrom TreeTools TipLabels
#' @param checks Logical specifying whether to perform basic sanity checks to
#' avoid crashes in C++.
#' @keywords internal
#' @seealso [`CalculateTreeDistance`]
#' @export
.TreeDistance <- function (Func, tree1, tree2, checks = TRUE, ...) {
  single1 <- inherits(tree1, 'phylo')
  labels1 <- TipLabels(tree1, single = TRUE)
  nTip <- length(labels1)
  
  single2 <- inherits(tree2, 'phylo')
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

.TreeDistanceOneMany <- function (Func, oneTree, manyTrees, tipLabels, 
                                  nTip = length(tipLabels), 
                                  FUN.VALUE = Func(manyTrees[[1]], oneTree,
                                                   nTip = nTip, 
                                                   tipLabels = tipLabels, ...),
                                  ...) {
  vapply(manyTrees, Func, oneTree, tipLabels = tipLabels,
                            nTip = nTip, ..., FUN.VALUE = FUN.VALUE)
}

.TreeDistanceManyMany <- function (Func, trees1, trees2, tipLabels, 
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

.CheckLabelsSame <- function (labelList) {
  nTip <- unique(vapply(labelList, length, 0L))
  if (length(nTip) != 1) {
    stop("All trees must contain the same number of leaves.")
  }
  tipLabel <- unique(lapply(labelList, sort))
  if (length(tipLabel) != 1L) {
    stop("All trees must bear identical labels. Found:\n >  ", 
         paste0(lapply(tipLabel, paste, collapse = ' '), 
               ' (', lapply(tipLabel, class), ')',
                collapse = '\n >  '))
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
Entropy <- function (...) {
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
#' library('TreeTools', quietly = TRUE, warn.conflicts = FALSE)
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
#' @importFrom stats dist
#' @export
CompareAll <- function (x, Func, FUN.VALUE = Func(x[[1]], x[[1]], ...),
                        ...) {
  nTree <- length(x)
  countUp <- seq_len(nTree - 1)
  i <- x[rep(countUp, rev(countUp))]
  j <- x[unlist(sapply(countUp, function (n) n + seq_len(nTree - n)))]
  
  ret <- vapply(seq_along(i), function (k) Func(i[[k]], j[[k]], ...), FUN.VALUE)
  
  .WrapReturn <- function (dists) {
    structure(dists,
              Size = nTree,
              Labels = names(x),
              Diag = FALSE,
              Upper = TRUE,
              class = 'dist')
  }
  
  # Return:
  if (length(FUN.VALUE) == 1) {
    .WrapReturn(ret)
  } else {
    structure(lapply(seq_len(dim(ret)[1]), 
                     function (i) .WrapReturn(unlist(ret[i, ]))),
              names = rownames(ret))
  }
}

#' Normalize information against total present in both starting trees
#' @param unnormalized Numeric value to be normalized.
#' @param tree1,tree2 Trees from which `unnormalized` was calculated
#' @param InfoInTree Function to calculate the information content of each tree
#' @param infoInBoth Numeric specifying information content of both trees
#' independently (optional)
#' @param how Method for normalization
#' @param Func Function that takes as inputs `tree1Info` and `tree2Info`, and
#' returns a normalizing constant against which to divide `unnormalized`.
#' @param \dots Additional parameters to `InfoInTree()` or `how`.
#' @keywords internal
#' @template MRS
#' @export
NormalizeInfo <- function (unnormalized, tree1, tree2, InfoInTree, 
                           infoInBoth = NULL,
                           how = TRUE, Combine = '+', ...) {
  
  CombineInfo <- function (tree1Info, tree2Info, Combiner = Combine) {
    if (length(tree1Info) == 1 || length(tree2Info) == 1) {
      mapply(Combiner, tree1Info, tree2Info)
    } else {
      ret <- outer(tree1Info, tree2Info, Combiner)
      if (inherits(unnormalized, 'dist')) ret[lower.tri(ret)] else ret
    }
  }
  
  if (is.logical(how)) {
    if (how == FALSE) {
      return (unnormalized)
    } else {
      if (is.null(infoInBoth)) {
        info1 <- InfoInTree(tree1, ...)
        info2 <- if (is.null(tree2)) info1 else InfoInTree(tree2, ...)
        infoInBoth <- CombineInfo(info1, info2)
      }
    }
  } else if (is.function(how)) {
    if (is.null(infoInBoth)) {
      info1 <- InfoInTree(tree1, ...)
      info2 <- if (is.null(tree2)) info1 else InfoInTree(tree2, ...)
      infoInBoth <- CombineInfo(info1, info2, Combiner = how)
    }
  } else {
    infoInBoth <- how
  }
  
  # Return:
  unnormalized / infoInBoth
}

#' List clades as text
#' @param splits,splits1,splits2 Logical matrices with columns specifying membership
#' of each corresponding matched clade.
#' @return `ReportMatching` returns a character vector describing each pairing 
#' in a matching.
#'   
#' @seealso [`VisualizeMatching`]
#' @template MRS
#' @keywords internal
#' @export
ReportMatching <- function (splits1, splits2, realMatch = TRUE) {
  paste(as.character(splits1), ifelse(realMatch, '=>', '..'),
        as.character(splits2))
}
