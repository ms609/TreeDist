#' @export
as.HPart_cpp <- function(tree, tipLabels) {
  UseMethod("as.HPart_cpp")
}

#' @export
as.HPart_cpp.HPart_cpp <- function(tree, tipLabels) {
  if (identical(tipLabels, TipLabels(tree))) {
    trees
  } else {
    stop("Relabelling not yet implemented")
  }
}

#' @export
as.HPart_cpp.list <- function(tree) {
  stop("Not yet implemented") # TODO
}

#' @export
as.HPart_cpp.phylo <- function(tree, tipLabels = TipLabels(tree)) {
  if (!identical(TipLabels(tree), tipLabels)) {
    tree <- RenumberTips(tree, tipLabels)
  }
  structure(build_hpart_from_phylo(tree), tip.label = tipLabels,
            class = "HPart_cpp")
}

#' @export
is.HPart_cpp <- function(x) {
  inherits(x, "HPart_cpp")
}

#' @export
NTip.HPart_cpp <- function(phy) {
  length(TipLabels(phy))
}

#' @export
print.HPart_cpp <- function(x, ...) {
  nTip <- NTip(x)
  tips <- TipLabels(x)
  cat("Hierarchical partition on", nTip, "leaves: ")
  if (nTip > 5) {
    cat(paste0(c(tips[1:2], "...", tips[length(tips) - 1:0]), collapse = ", "))
  } else {
    cat(paste0(tips, collapse = ", "))
  }
}

#' @export
as.phylo.HPart_cpp <- function(x, ...) {
  edge <- hpart_to_edge(x)
  labels <- TipLabels(x)
  nNode <- dim(edge)[[1]] - length(labels) + 1
  structure(list(edge = edge, Nnode = nNode, tip.label = labels),
            class = "phylo",
            order = "cladewise")
}

#' @export
plot.HPart_cpp <- function(x, ...) {
  plot(as.phylo(x), ...)
}

#' @export
clone <- function(x, ...) UseMethod("clone")

#' @export
clone.HPart_cpp <- function(x, tipLabel = attr(x, "tip.label")) {
  structure(clone_hpart(x), tip.label = tipLabel,
            class = "HPart_cpp")
}

#' @importFrom TreeTools MatchStrings
#' @export
RenumberTips.HPart_cpp <- function(tree, tipOrder) {
  startOrder <- TipLabels(tree)
  newOrder <- MatchStrings(TipLabels(tipOrder, single = TRUE), startOrder)
  
  if (!identical(newOrder, startOrder)) {
    newIndices <- match(newOrder, startOrder)
    if (any(is.na(newIndices))) {
      stop("Tree labels ", paste0(startOrder[is.na(newIndices)], collapse = ", "),
           " missing from `tipOrder`")
    }
    tree <- clone(tree, newOrder)
    relabel_hpart(tree, newIndices)
    # Return:
    tree
  } else {
    clone(tree)
  }
}