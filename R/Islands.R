#' Find islands from distance matrix
#' 
#' `Islands()` assigns a set of objects to islands, such that all elements
#' within an island can form a connected graph in which each edge is no longer
#' than `threshold` distance units \insertRef{Silva2021}{TreeDist}.
#' 
#' @inheritParams SpectralEigens
#' @param threshold Elements greater than `threshold` distance units will not be
#' assigned to the same island.
#' @param dense Logical; if `FALSE`, each island will be named according to the
#' index of its lowest-indexed member; if `TRUE`, each island will be numbered
#' sequentially from `1`, ordered by the index of the lowest-indexed member.
#' @param smallest Integer; Islands comprising no more than `smallest` elements
#' will be assigned to island `NA`.
#' @return `Islands()` returns a vector listing the island to which
#' each element is assigned.
#' @references \insertAllCited{}
#' @examples
#' library("TreeTools", quietly = TRUE)
#' # Generate a set of trees
#' trees <- as.phylo(as.TreeNumber(BalancedTree(16)) + c(-(40:20), 70:105), 16)
#' 
#' # Calculate distances between trees
#' distances <- ClusteringInfoDist(trees)
#' summary(distances)
#' 
#' # Assign trees to islands
#' isle <- Islands(distances, quantile(distances, 0.1))
#' table(isle)
#' 
#' # Indicate island membership on 2D mapping of tree distances
#' mapping <- cmdscale(distances, 2)
#' plot(mapping, col = isle + 1,
#'      asp = 1, # Preserve aspect ratio - do not distort distances
#'      ann = FALSE, axes = FALSE, # Don't label axes: dimensions are meaningless)
#'      pch = 16 # Plotting character: Filled circle
#' )
#' 
#' # Compare strict consensus with island consensus trees
#' oPar <- par(mfrow = c(2, 2), mai = rep(0.1, 4))
#' plot(Consensus(trees), main = "Strict")
#' plot(Consensus(trees[isle == 1]), edge.col = 2, main = "Island 1")
#' plot(Consensus(trees[isle == 2]), edge.col = 3, main = "Island 2")
#' plot(Consensus(trees[isle == 3]), edge.col = 4, main = "Island 3")
#' 
#' # Restore graphical parameters
#' par(oPar)
#' @template MRS
#' @family tree space functions
#' @export
Islands <- function(D, threshold, dense = TRUE, smallest = 0) {
  linked <- as.matrix(D) <= threshold
  n <- dim(linked)[[1]]
  ret <- integer(n)
  i <- 1
  repeat {
    links <- seq_len(n) == i
    repeat {
      nowLinked <- colSums(linked[links, , drop = FALSE]) > 0
      if (any(nowLinked[!links])) {
        # Added to island
        links <- nowLinked
      } else {
        break
      }
    }
    ret[links] <- i
    i <- which.min(ret)
    if (ret[[i]]) break
  }
  tab <- table(ret, dnn = NULL)
  ret[ret %in% as.integer(names(tab)[tab < smallest])] <- NA
  
  if (dense) {
    as.integer(factor(rank(ret, ties.method = "min", na.last = "keep")))
  } else {
    ret
  }
}