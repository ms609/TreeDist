
#' @export
Fitch <- function (tree, data, at = NULL) {
  if (is.null(at)) at <- attributes(data)
  n.char  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  info <- at$info.amounts
  if (is.null(at$order) || at$order == "cladewise") tree <- reorder(tree, "postorder")
  tree.edge <- tree$edge
  parent <- tree.edge[, 1]
  child <- tree.edge[, 2]
  tip.label <- tree$tip.label
  n.edge <- length(parent)
  max.node <- parent[1] #max(parent)
  n.tip <- length(tip.label)
  n.node <- max.node - n.tip
  inapp <- at$inapp.level
  parent.of <- parent[match((n.tip + 2L):max.node, child )]
  allNodes <- (n.tip + 1L):max.node
  child.of <- child [c(match(allNodes, parent),
                       length(parent) + 1L - match(allNodes, rev(parent)))]
  fitch <- .Call("FITCH", data[, tree$tip.label], as.integer(n.char),
        as.integer(parent), as.integer(child), as.integer(n.edge),
        as.double(weight), as.integer(max.node), as.integer(n.tip), PACKAGE='phangorn')
  return(fitch[[2]])
#
#  TODO: link library "inapplicable" to support inapplicable data...
}