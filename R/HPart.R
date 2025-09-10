#' Hierarchical partition structure
#' 
#' @name HPart
#' @export
as.HPart <- function(tree, tipLabels) {
  UseMethod("as.HPart")
}

#' @export
#' @rdname HPart
as.HPart.HPart <- function(tree, tipLabels = NULL) {
  if (is.null(tipLabels) || identical(tipLabels, TipLabels(tree))) {
    tree
  } else {
    RenumberTips(tree, tipLabels)
  }
}

#' @param tree hierarchical list-of-lists (leaves = integers 1..n)
#' @export
as.HPart.list <- function(tree, tipLabels = NULL) {
  # Flatten to verify leaves
  leaves <- unlist(tree, recursive = TRUE)
  if (!all(is.numeric(leaves)) || any(leaves != as.integer(leaves))) {
    stop("All leaves must be integers.")
  }
  tree <- rapply(tree, as.integer, how = "replace")
  if (0 %in% leaves) {
    tree <- rapply(tree, function(x) x + 1L, how = "replace")
    leaves <- leaves + 1
  }
  n_tip <- length(unique(leaves))
  expected <- seq_len(n_tip)
  if (!setequal(leaves, expected)) {
    stop("Leaves must contain all integers 1..n without gaps")
  }
  
  hpart_ptr <- build_hpart_from_list(tree, n_tip)
  ret <- structure(hpart_ptr, tip.label = as.character(expected), class = "HPart")
  if (!is.null(tipLabels)) {
    RenumberTips(ret, tipLabels)
  }
  ret
}

#' @export
as.HPart.phylo <- function(tree, tipLabels = TipLabels(tree)) {
  if (!identical(TipLabels(tree), tipLabels)) {
    tree <- RenumberTips(tree, tipLabels)
  }
  structure(build_hpart_from_phylo(tree), tip.label = tipLabels,
            class = "HPart")
}

#' @export
is.HPart <- function(x) {
  inherits(x, "HPart")
}

#' @export
NTip.HPart <- function(phy) {
  length(TipLabels(phy))
}

#' @export
print.HPart <- function(x, ...) {
  nTip <- NTip(x)
  tips <- TipLabels(x)
  cat("Hierarchical partition on", nTip, "leaves: ")
  if (nTip > 5) {
    cat(paste0(c(tips[1:2], "...", tips[length(tips) - 1:0]), collapse = ", "))
  } else {
    cat(paste0(tips, collapse = ", "))
  }
}

#' @importFrom ape as.phylo
#' @export
as.phylo.HPart <- function(x, ...) {
  edge <- hpart_to_edge(x)
  labels <- TipLabels(x)
  nNode <- dim(edge)[[1]] - length(labels) + 1
  structure(list(edge = edge, Nnode = nNode, tip.label = labels),
            class = "phylo",
            order = "cladewise")
}

#' @export
plot.HPart <- function(x, ...) {
  plot(as.phylo(x), ...)
}

#' @export
clone <- function(x, ...) UseMethod("clone")

#' @export
clone.HPart <- function(x, tipLabel = attr(x, "tip.label")) {
  structure(clone_hpart(x), tip.label = tipLabel,
            class = "HPart")
}

#' @importFrom TreeTools MatchStrings
#' @export
RenumberTips.HPart <- function(tree, tipOrder) {
  startOrder <- TipLabels(tree)
  newOrder <- MatchStrings(TipLabels(tipOrder, single = TRUE), startOrder)
  
  if (!identical(newOrder, startOrder)) {
    newIndices <- match(newOrder, startOrder)
    if (any(is.na(newIndices))) {
      stop("Tree labels ", paste0(startOrder[is.na(newIndices)], collapse = ", "),
           " missing from `tipOrder`")
    }
    tree <- clone(tree, newOrder)
    relabel_hpart(tree, newIndices - 1L)
    # Return:
    tree
  } else {
    clone(tree)
  }
}
