#' Kendall&ndash;Colijn distance
#' 
#' Calculate the Kendall&ndash;Colijn tree distance, a measure related to the 
#' path difference. 
#' 
#' The Kendall&ndash;Colijn distance works by measuring, for each pair of leaves,
#' the distance from the most recent
#' common ancestor of those leaves and the root node.  For a given tree, this 
#' produces a vector of values recording the distance-from-the-root of each
#' most recent common ancestor of each pair of leaves.
#' 
#' Two trees are compared by taking the Euclidian distance between the
#' respective vectors.  This is calculated by taking the square root of the sum 
#' of the squares of the differences between the vectors.
#' 
#' This metric emphasizes the position of the root; the path difference 
#' instead measures the distance of the last common ancestor of each pair
#' of leaves from the leaves themselves, i.e. the length of the path from one 
#' leaf to another.
#' 
#' @template tree12ListParams
#' @param Vector Function convering a tree to a numeric vector. 
#' `KCVector`, the default, returns the number of edges between the common
#' ancestor of each pair of leaves and the root of the tree (per
#' Kendall & Colijn 2016).
#' `PathVector` returns the number of edges between each pair of leaves (per
#' Steel & Penny 1993).
#' `SplitVector` returns the size of the smallest split that contains each
#' pair of leaves (per Smith, forthcoming).
#' 
#' @templateVar returns `KendallColijn()` returns
#' @template distReturn
#' 
#' @examples 
#' KendallColijn(TreeTools::BalancedTree(8), TreeTools::PectinateTree(8))
#'
#' set.seed(0)
#' KendallColijn(TreeTools::BalancedTree(8), lapply(rep(8, 3), ape::rtree))
#' KendallColijn(lapply(rep(8, 4), ape::rtree))
#' 
#' KendallColijn(lapply(rep(8, 4), ape::rtree), Vector = SplitVector)
#' 
#' # Notice that changing tree shape close to the root results in much
#' # larger differences
#' tree1 <- ape::read.tree(text = "(a, (b, (c, (d, (e, (f, (g, h)))))));")
#' tree2 <- ape::read.tree(text = "(a, ((b, c), (d, (e, (f, (g, h))))));")
#' tree3 <- ape::read.tree(text = "(a, (b, (c, (d, (e, ((f, g), h))))));")
#' trees <- c(tree1, tree2, tree3)
#' KendallColijn(trees)
#' KendallColijn(trees, Vector = SplitVector)
#' @template MRS
#' 
#' @seealso [`treespace::treeDist`](https://CRAN.R-project.org/package=treespace/vignettes/introduction.html)
#' is a more sophisticated, if more cumbersome, implementation that supports 
#' lambda > 0, i.e. use of edge lengths in tree comparison.
#' 
#' @references \insertRef{Kendall2016}{TreeDist}
#' 
#' @references \insertRef{SmithSpace}{TreeDist}
#' @family tree distances
#' @encoding UTF-8
#' @export
KendallColijn <- function (tree1, tree2 = NULL, Vector = KCVector) {
  
  FunValue <- function (nTip) double(nTip * (nTip - 1L) / 2L)
  
  EuclidianDistance <- function (x) sqrt(sum(x * x))
  
  if (inherits(tree1, 'phylo')) {
    if (inherits(tree2, 'phylo')) {
      if (length(tree1$tip.label) != length(tree2$tip.label) || 
          length(setdiff(tree1$tip.label, tree2$tip.label)) > 0) {
        stop("Leaves must bear identical labels.")
      }
      EuclidianDistance(Vector(tree1) - Vector(tree2))
    } else {
      if (is.null(tree2)) {
        0
      } else {
        apply(Vector(tree1) - vapply(tree2, Vector,
                                     FunValue(length(tree1$tip.label))),
              2L, EuclidianDistance)
      }
    }
  } else {
    if (inherits(tree2, 'phylo')) {
      apply(Vector(tree2) - vapply(tree1, Vector,
                                   FunValue(length(tree2$tip.label))),
            2L, EuclidianDistance)
    } else if (is.null(tree2)) {
      
      treeVec <- vapply(tree1, Vector, FunValue(length(tree1[[1]]$tip.label)))
      nTree <- length(tree1)
      ret <- matrix(0, nTree, nTree)
      is <- matrix(seq_len(nTree), nTree, nTree)
      
      ret <- structure(class = 'dist', Size = nTree,
                       diag = FALSE, upper = FALSE,
                       apply(treeVec[, t(is)[lower.tri(is)]] - 
                               treeVec[, is[lower.tri(is)]], 2,
                             EuclidianDistance))
      # Return:
      ret
    }
      else {
      vector1 <- vapply(tree1, Vector, FunValue(length(tree1[[1]]$tip.label)))
      vector2 <- vapply(tree2, Vector, FunValue(length(tree2[[1]]$tip.label)))
      apply(vector2, 2, function(i) 
        apply(vector1, 2, function (j) 
          EuclidianDistance(i - j)))
    }
  }
}

#' @describeIn KendallColijn Creates a vector that characterises a rooted tree,
#' as described in Kendall & Colijn (2016).
#' @param tree A tree of class \code{\link[ape:read.tree]{phylo}}.
#' @importFrom TreeTools AllAncestors Preorder
#' @export
KCVector <- function (tree) {
  tree <- Preorder(tree)
  edge <- tree$edge
  parent <- edge[, 1L]
  child <- edge[, 2L]
  root <- parent[1]
  nTip <- root - 1L
  tipOrder <- order(tree$tip.label)
  
  tipI <- rep(tipOrder, nTip - seq_len(nTip))
  tipJ <- matrix(tipOrder, nrow = nTip, ncol = nTip)
  tipJ <- tipJ[lower.tri(tipJ)]
  
  ancestors <- AllAncestors(parent, child)
  
  mrca <- mapply(function(anc1, anc2) max(intersect(anc1, anc2)), 
                 ancestors[tipI], ancestors[tipJ])
  
  rootDist <- vapply(ancestors, length, integer(1))
  structure(rootDist[mrca], Size = nTip, class = 'dist')
}

#' @describeIn KendallColijn Creates a vector reporting the path length between
#' each pair of leaves, per the path metric of Steel & Penny (1993).
#' @importFrom TreeTools AllAncestors Preorder
#' @export
PathVector <- function (tree) {
  tree <- Preorder(tree)
  edge <- tree$edge
  parent <- edge[, 1L]
  child <- edge[, 2L]
  root <- parent[1]
  nTip <- root - 1L
  tipAncs <- seq_len(nTip)
  tipOrder <- order(tree$tip.label)
  
  tipI <- rep(tipOrder, nTip - seq_len(nTip))
  tipJ <- matrix(tipOrder, nrow = nTip, ncol = nTip)
  tipJ <- tipJ[lower.tri(tipJ)]
  
  ancestors <- AllAncestors(parent, child)
  
  pathLength <- mapply(function(anc1, anc2) {
    mrca <- max(intersect(anc1, anc2))
    sum(anc1 >= mrca, anc2 >= mrca,
        -(mrca == root) # don't count root edge twice
        )
  }, ancestors[tipI], ancestors[tipJ])
  
  structure(pathLength, Size = nTip, class = 'dist')
}

#' @describeIn KendallColijn Creates a vector reporting the smallest split
#' containing each pair of leaves, per the metric proposed in Smith
#' (forthcoming).
#' @importFrom TreeTools as.Splits
#' @export
SplitVector <- function (tree) {
  tipLabel <- tree$tip.label
  nTip <- length(tipLabel)
  splits <- as.logical(as.Splits(tree, tipLabel[order(tipLabel)]))
  splits <- rbind(splits, !splits)
  inSplit <- rowSums(splits)
  
  tipI <- rep(seq_len(nTip), nTip - seq_len(nTip))
  tipJ <- matrix(seq_len(nTip), nrow = nTip, ncol = nTip)
  tipJ <- tipJ[lower.tri(tipJ)]
  
  smallestSplit <- mapply(function(tip1, tip2) {
    min(inSplit[splits[, tip1] & splits[, tip2]], nTip - 1L)
  }, tipI, tipJ)
  
  structure(smallestSplit, Size = nTip, class = 'dist')
}
