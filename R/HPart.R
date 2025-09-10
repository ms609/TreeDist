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
