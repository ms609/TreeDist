
#' @keywords internal
#' @export
Min <- function (x, inappLevel) {
  if (length(inappLevel)) return(sum(2^(c(0:(inappLevel-2), inappLevel:12)) %in% unique(x)))
  return (sum(2^(0:12) %in% unique(x)))
}


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
#' @param nEdge number of edges (calcluated from length(parent) if not supplied)
#' @return a logical vector stating whether each edge in turn is a descendant of the speficied edge
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

#' @keywords internal
#' @export
AncestorEdge <- function (edge, parent, child) child == parent[edge]

#' Descendant Edges
#'
#' Quickly identifies edges that are 'ancestral' to a particular edge in a tree
#'
#' @param edge number of the edge whose child edges are required
#' @template treeParent
#' @template treeChild
#' @param stopAt number of the edge at which the search should terminate; defaults to the root edges
#' @param nEdge number of edges (calcluated from length(parent) if not supplied)
#' @return a logical vector stating whether each edge in turn is a descendant of the speficied edge
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
#' 
#' @author Martin R. Smith 
#' @importFrom ape rtree
#' @export
RandomTree <- function (dataset) rtree(length(dataset), tip.label=names(dataset), br=NULL)
