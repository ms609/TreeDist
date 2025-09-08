#' @export
as.HPart <- function(tree) {
  UseMethod("as.HPart")
}

#' @export
as.HPart.HPart <- function(tree) tree

#' @export
as.HPart.list <- function(tree) {
  structure(tree, class = "HPart")
}

#' @export
as.HPart.phylo <- function(tree) {
  # Ensure tree is rooted and binary (ape usually handles this)
  edge <- Preorder(tree$edge)
  tips <- tree$tip.label
  nTip <- length(tips)
  
  # Build adjacency list
  children <- vector("list", nTip + tree$Nnode)
  for (i in seq_len(nrow(edge))) {
    parent <- edge[i, 1]
    child  <- edge[i, 2]
    children[[parent]] <- c(children[[parent]], child)
  }
  
  # Recursive builder
  .Build <- function(node) {
    kids <- children[[node]]
    if (length(kids) == 0) {
      list(tips[[node]])
    } else {
      leaves <- kids <= nTip
      if (all(leaves)) {
        as.list(tips[kids])
      } else {
        lapply(children[[node]], .Build)
      }
    }
  }
  
  root <- nTip + 1
  structure(.Build(root), class = "HPart")
}

.ValidPartition <- function(x) {
  if (all(vapply(x, is.list, logical(1)))) {
    all(vapply(x, .ValidPartition, logical(1)))
  } else {
    all(vapply(x, is.character, logical(1))) ||
      all(vapply(x, is.numeric, logical(1)))
  }
}

# Replicates check(hp)
#' @source https://github.com/jipphysics/hit/blob/master/hit.ipynb
#' @export
is.HPart <- function(x) {
  inherits(x, "HPart") && .ValidPartition(x)
}
