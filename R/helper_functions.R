#' Quick sample
#' 
#' Faster than inbuilt sample because it avoids some checks
#' @param x vector to sample
#' @param len length of vector
#' @keywords internal
#' @export
SampleOne <- function (x, len = length(x)) x[sample.int(len, 1L, FALSE, NULL, FALSE)]

Assert <- function (statement) if (!statement) stop(deparse(statement), " is FALSE")

#' Descendant Edges
#'
#' Quickly identifies edges that are 'descended' from a particular edge in a tree
#'
#' @param edge number of the edge whose child edges are required
#' @template treeParent
#' @template treeChild
#' @param nEdge number of edges (calculated from length(parent) if not supplied)
#' @return a logical vector stating whether each edge in turn is a descendant of the specified edge
#'         (or the edge itself)
#' @export
DescendantEdges <- function (edge, parent, child, nEdge = length(parent)) {
  ret <- logical(nEdge)
  edgeSister <- match(parent[edge], parent[-edge])
  if (edgeSister >= edge) {
    # edgeSister is really 1 higher than you think, because we knocked out edge 'edge' in the match
    ret[edge:edgeSister] <- TRUE
    return(ret)
  } else {
    nextEdge <- edge
    revParent <- rev(parent)
    repeat {
      if (revDescendant <- match(child[nextEdge], revParent, nomatch=FALSE)) {
        nextEdge <- 1 + nEdge - revDescendant
      } else break;
    }
    ret[edge:nextEdge] <- TRUE 
    return(ret)
  }
}

#' Ancestral edge
#'
#' @param edge Number of an edge
#' @param parent parent nodes (given by phylo.object$edge[, 1])
#' @param child  child nodes (given by phylo.object$edge[, 2])
#' @return a logical vector identifying whether each edge is the edge that is ancestral to the given edge.
#' @keywords internal
#' @export
AncestorEdge <- function (edge, parent, child) child == parent[edge]

#' EdgeAncestry
#'
#' Descendant Edges
#'
#' Quickly identifies edges that are 'ancestral' to a particular edge in a tree
#'
#' @param edge number of the edge whose child edges are required
#' @template treeParent
#' @template treeChild
#' @param stopAt number of the edge at which the search should terminate; defaults to the root edges
#' @return a logical vector stating whether each edge in turn is a descendant of the specified edge
#'
#' @author Martin R. Smith
#' @export
EdgeAncestry <- function (edge, parent, child, stopAt = (parent==min(parent))) {
  ret <- edge <- AncestorEdge(edge, parent, child)
  repeat {
    if (any(ret[stopAt])) return(ret)
    ret[edge <- AncestorEdge(edge, parent, child)] <- TRUE    
  }
}

#' Generate random tree topology from dataset
#' 
#' @param dataset A dataset in \code{\link[phangorn]{phyDat}} format
#' @param root Taxon to use as root (if desired; FALSE otherwise)
#' 
#' @author Martin R. Smith 
#' @importFrom ape rtree
#' @importFrom ape root
#' @export
RandomTree <- function (dataset, root = FALSE) {
  tree <- rtree(length(dataset), tip.label=names(dataset), br=NULL)
  return (if (root != FALSE) root(tree, root, resolve.root=TRUE) else tree)
}

#' Force taxa to form an outgroup
#'
#' Given a tree or a list of taxa, rearrange the ingroup and outgroup taxa such that the two
#' are sister taxa across the root, without altering the relationships within the ingroup
#' or within the outgroup.
#'
#' @param tree either a tree of class \code{phylo}, or a character vector listing the names of 
#'        all the taxa in the tree, from which a random tree will be generated.
#' @param outgroup a vector containing the names of taxa to include in the outgroup
#'
#' @return a tree where all outgroup taxa are sister to all remaining taxa, 
#'         otherwise retaining the topology of the ingroup.
#' @author Martin R. Smith
#' @importFrom ape rtree
#' @importFrom ape root drop.tip bind.tree
#' @export
EnforceOutgroup <- function (tree, outgroup) {
  if (class(tree) == 'phylo') {
    taxa <- tree$tip.label
  } else if (class(tree) == 'character') {    
    tree <- root(rtree(length(taxa), tip.label=taxa, br=NULL), taxa[1], resolve.root=TRUE)
  } else {
    stop ("tree must be of class phylo")
  }
  
  if (length(outgroup) == 1) return (root(tree, outgroup, resolve.root=TRUE))
  
  ingroup <- taxa[!(taxa %in% outgroup)]
  if (!all(outgroup %in% taxa) || length(ingroup) + length(outgroup) != length(taxa)) {
    stop ("All outgroup taxa must occur in speficied taxa")
  }
  
  ingroup.branch <- drop.tip(tree, outgroup)
  outgroup.branch <- drop.tip(tree, ingroup)
  
  result <- root(bind.tree(outgroup.branch, ingroup.branch, 0, 1), outgroup, resolve.root=TRUE)
  RenumberTips(Renumber(result), taxa)
}
