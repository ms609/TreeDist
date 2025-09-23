#' Hierarchical partition structure
#' 
#' A structure of class `HPart` comprises a pointer to a C++ representation of
#' hierarchical partitions, with the attribute `tip.label` recording the
#' character labels of its leaves.  `HPart` objects with identical tip labels
#' can be compared using [`HierarchicalMutualInfo()`].
#' 
#' 
#' An `HPart` object may be created from various representations of hierarchical
#' structures:
#' 
#' - a tree of class `phylo`
#' - A hierarchical list of lists, in which elements are represented by integers
#'   1\dots{}n
#' - A vector, which will be interpreted as a flat structure
#'  in which all elements bearing the same label are assigned to the same cluster
#' 
#' @param tree An object to convert to an HPart structure, in a supported format
#' (see details)
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

#' @rdname HPart
#' @export
as.HPart.default <- function(tree, tipLabels = NULL) {
  if (is.null(dim(tree))) {
    structure(build_hpart_from_list(
      lapply(unique(tree), function(x) as.list(which(tree == x))),
      length(tree)),
      tip.label = seq_along(tree),
      class = "HPart"
    )
  } else {
    stop("no applicable method for 'as.HPart' applied to an object of class \"",
         paste(class(tree), collapse = "\", \""), "\"")
  }
}


#' @rdname HPart
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
  if (!isTRUE(all.equal(sort(leaves), expected))) {
    stop("Leaves must contain each integer 1..n exactly once")
  }
  
  hpart_ptr <- build_hpart_from_list(tree, n_tip)
  ret <- structure(hpart_ptr, tip.label = as.character(expected), class = "HPart")
  if (!is.null(tipLabels) && !is.list(tipLabels)) {
    RenumberTips(ret, tipLabels)
  }
  ret
}

#' @export
#' @inheritParams TreeTools::as.ClusterTable
#' @rdname HPart
as.HPart.phylo <- function(tree, tipLabels = TipLabels(tree)) {
  if (!identical(TipLabels(tree), tipLabels)) {
    tree <- RenumberTips(tree, tipLabels)
  }
  structure(build_hpart_from_phylo(tree), tip.label = tipLabels,
            class = "HPart")
}

#' @rdname HPart
#' @export
is.HPart <- function(x) {
  inherits(x, "HPart")
}

#' @export
NTip.HPart <- function(phy) {
  length(TipLabels(phy))
}

#' @rdname HPart
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

#' @rdname HPart
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

#' @rdname HPart
#' @param x `HPart` object to plot
#' @param \dots Additional arguments to \code{\link[ape:plot.phylo]{plot.phylo}}
#' @export
plot.HPart <- function(x, ...) {
  plot(as.phylo(x), ...)
}

#' Clone / duplicate an object
#' `clone()` physically duplicates objects
#' @param x the object to be cloned
#' @param \dots additional parameters for methods
#' @return `clone()` typically returns an object of the same class and "value"
#' as the input `x`.
#' @export
clone <- function(x, ...) UseMethod("clone")

#' @template MRS
#' @rdname clone
#' @inheritParams TreeTools::as.ClusterTable
#' @export
clone.HPart <- function(x, tipLabels = attr(x, "tip.label"), ...) {
  structure(clone_hpart(x), tip.label = tipLabels,
            class = "HPart")
}

#' @importFrom TreeTools MatchStrings
#' @inheritParams TreeTools::RenumberTips
#' @export
RenumberTips.HPart <- function(tree, tipOrder) {
  startOrder <- TipLabels(tree)
  newOrder <- MatchStrings(TipLabels(tipOrder, single = TRUE), startOrder)
  
  if (!identical(newOrder, startOrder)) {
    if (length(newOrder) != length(startOrder)) {
      stop("Tree labels ", paste0(setdiff(startOrder, tipOrder), collapse = ", "),
           " missing from `tipOrder`")
    }
    newIndices <- match(newOrder, startOrder)
    tree <- clone(tree, newOrder)
    relabel_hpart(tree, newIndices - 1L)
    # Return:
    tree
  } else {
    clone(tree)
  }
}
