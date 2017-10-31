#' Renumber a tree's nodes and tips
#'
#' \code{Renumber} numbers the nodes and tips in a tree to conform with the phylo standards.
#'
#' @template treeParam 
#' 
#' @examples
#' library('ape')
#' tree <- rtree(10)
#' Renumber (tree)
#' 
#' @return This function returns a tree of class \code{phylo}
#'   
#' @author Martin R. Smith
#' @export
Renumber <- function (tree) {
  tree   <- Postorder(tree)
  edge   <- tree$edge
  nTip   <- length(tree$tip.label)
  parent <- edge[, 1L]
  child  <- edge[, 2L]
  NODES  <- child > nTip
  TIPS   <- !NODES
  nNode  <- sum(NODES) + 1 # Root node has no edge leading to it, so add 1
  
  tip <- child[TIPS]
  name <- vector("character", length(tip))
  name[1:nTip] <- tree$tip.label[tip]
  tree$tip.label <- name
  child[TIPS] <- 1:nTip
  
  old.node.number <- unique(parent)
  new.node.number <- rev(nTip + seq_along(old.node.number))
  renumbering.schema <- integer(nNode)
  renumbering.schema[old.node.number - nTip] <- new.node.number
  child[NODES] <- renumbering.schema[child[NODES] - nTip]
  nodeseq <- (1L:nNode) * 2L
  parent <- renumbering.schema[parent - nTip]
  
  tree$edge[,1] <- parent
  tree$edge[,2] <- child
  Cladewise(tree)
}

#' SingleTaxonTree
#'
#' Single taxon tree
#'
#' Create a phylogenetic 'tree' that comprises a single taxon.
#'
#' @usage SingleTaxonTree(label)
#' @param  label a character vector specifying the label of the tip.
#' @return This function returns a \code{phylo} object containing a single tip with the specified label.
#' @examples SingleTaxonTree('Homo_sapiens')
#' @keywords  tree 
#' @export
SingleTaxonTree <- function (label) {
  res <- list(edge=matrix(c(2L,1L), 1, 2), tip.label=label, Nnode=1L)
  class(res) <- 'phylo'
  res
}

#' Extract subtree
#'
#' @description Safely extracts a clade from a phylogenetic tree.
#' @usage Subtree(tree, node)
#' 
#' 
#' @template preorderTreeParam
#' @param node The number of the node at the base of the clade to be extracted.
#' 
#' @details
#' Modified from the \pkg{ape} function \code{\link{extract.clade}}, which sometimes behaves erratically.  
#' Unlike extract.clade, this function supports the extraction of 'clades' that constitute a single tip.
#' 
#' @return This function returns a tree of class \code{phylo} that represents a clade 
#'         extracted from the original tree.
#'
#' @examples{
#' tree <- Preorder(ape::rtree(20, br=NULL))
#' plot(tree); ape::nodelabels(); ape::nodelabels(33, 33, bg='yellow'); dev.new()
#' plot(Subtree(tree, 33))
#' }
#' 
#' @author Martin R. Smith
#' @export
Subtree <- function (tree, node) {
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'preorder') stop("Tree must be in preorder")
  tipLabel <- tree$tip.label
  nTip <- length(tipLabel)
  if (node <= nTip) return(SingleTaxonTree(tipLabel[node]))
  if (node == nTip + 1L) return(tree)

  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  subtreeParentEdge <- match(node, child)
  keepEdge <- DescendantEdges(subtreeParentEdge, parent, child)
  keepEdge[subtreeParentEdge] <- FALSE
  
  edge <- edge[keepEdge, ]
  edge1 <- edge[, 1]
  edge2 <- edge[, 2]
  
  isTip <- edge2 <= nTip
  tips  <- edge2[isTip]
  new.nTip <- length(tips)
  name <- character(new.nTip)
  name[order(tips)] <- tipLabel[tips]
  edge2[isTip] <- order(tips)
    
  ## renumber nodes:
  nodeAdjust <- new.nTip + 1 - node
  edge2[!isTip] <- edge2[!isTip] + nodeAdjust
  edge[, 1] <- edge1 + nodeAdjust
  edge[, 2] <- edge2
  newtree <- list(
    tip.label = name,
    Nnode = dim(edge)[1] - new.nTip + 1L,
    edge = edge
  )
  class(newtree)<-'phylo'
  newtree
}

#' Add a tip to a phylogenetic tree
#' 
#' \code{AddTip} adds a tip to a phylogenetic tree at a specified location.
#'
#' \code{AddTip} extends \code{\link{bind.tree}}, which cannot handle single-taxon trees.
#'
#' @usage AddTip(tree, where, label)
#' 
#' @template treeParam 
#' @param where The node or tip that should form the sister taxon to the new node.  To add a new tip at the root, use "where = 0";
#' @param label A character string providing the label the new tip.
#' 
#' @return This function returns a tree of class \code{phylo} with an additional tip at the desired location.
#' 
#' @author Martin R. Smith
#' 
#' @seealso \code{\link{bind.tree}}
#' @seealso \code{\link{nodelabels}}
#' 
#' @examples {
#'   library('ape')
#'   plot(tree <- rtree(10, br=NULL)); nodelabels(); nodelabels(15, 15, bg='green'); dev.new()
#'   plot(AddTip(tree, 15, 'NEW_TIP'))
#' }
#' @keywords tree 
#' 
#' @export
AddTip <- function (tree, where, label) {
  nTip <- length(tree$tip.label)
  nNode <- tree$Nnode
  ROOT <- nTip + 1L
  if (where < 1L) where <- ROOT
  new.tip.number <- nTip + 1L
  tree.edge <- tree$edge
  
  ## find the row of 'where' before renumbering
  if (where == ROOT) case <- 1 else {
      insertion.edge <- which(tree.edge[, 2] == where)
      case <- if (where <= nTip) 2 else 3
  }
  ## case = 1 -> y is bound on the root of x
  ## case = 2 -> y is bound on a tip of x
  ## case = 3 -> y is bound on a node of x

### because in all situations internal nodes need to be
### renumbered, they are changed to negatives first, and
### nodes eventually added will be numbered sequentially
  nodes <- tree.edge > nTip
  tree.edge[nodes] <- -(tree.edge[nodes] - nTip)  # -1, ..., -nTip
  next.node <- -nNode - 1L
  ROOT <- -1L # This may change later
  
  switch(case, { # case = 1 -> y is bound on the root of x
      tree.edge <- rbind(c(next.node, tree.edge[1]), tree.edge, c(next.node, new.tip.number))
      ROOT <- next.node
    }, { # case = 2 -> y is bound on a tip of x
      tree.edge[insertion.edge, 2] <- next.node
      tree.edge <- rbind(tree.edge[1:insertion.edge, ], c(next.node, where), c(next.node, new.tip.number), tree.edge[-(1:insertion.edge), ])
    }, { # case = 3 -> y is bound on a node of x
      tree.edge <- rbind(tree.edge[1:insertion.edge, ], c(next.node, tree.edge[insertion.edge, 2]), tree.edge[-(1:insertion.edge), ])
      tree.edge[insertion.edge, 2] <- next.node
      insertion.edge <- insertion.edge + 1L
      tree.edge <- rbind(tree.edge[1:insertion.edge, ], c(next.node, new.tip.number), tree.edge[-(1:insertion.edge), ])
    }
  )
  tree$tip.label <- c(tree$tip.label, label)
  tree$Nnode <- nNode <- nNode + 1L
  
  ## renumber nodes:
  new.numbering <- integer(nNode)
  new.numbering[-ROOT] <- new.tip.number + 1L
  second.col.nodes <- tree.edge[, 2] < 0
  ## executed from right to left, so newNb is modified before x$edge:
  tree.edge[second.col.nodes, 2] <- new.numbering[-tree.edge[second.col.nodes, 2]] <- new.tip.number + 2:nNode
  tree.edge[, 1] <- new.numbering[-tree.edge[, 1]]

  tree$edge <- tree.edge
  tree
  
}

#' List all ancestral nodes
#'
#' \code{AllAncestors} lists ancestors of each parent node in a tree
#'
#' Note that the tree's edges must be listed in an order whereby each entry in tr$edge[, 1] (with
#' the exception of the root) has appeared already in tr$edge[, 2]
#'
#' @template treeParent
#' @template treeChild
#' 
#' @examples
#' 
#' tr <- ape::rtree(20, br=NULL)
#' edge <- tr$edge
#' AllAncestors(edge[, 1], edge[, 2])
#' 
#' @return This function returns a list. Entry i contains a vector containing, in order,
#' the nodes encountered when traversing the tree from node i to the root node.  The last 
#' entry of each member of the list will therefore be the root node, with the exeption of the 
#' entry for the root node itself, which will be NULL.
#'   
#' @author Martin R. Smith
#' @export
AllAncestors <- function (parent, child) {
  res <- vector("list", max(parent))
  for (i in seq_along(parent)) {
    pa <- parent[i]
    res[[child[i]]] <- c(pa, res[[pa]])
  }
  res
}

#' Clade sizes
#' @template treeParam
#' @param nodes whose descendants should be returned
#' @return the number of nodes (including tips) that are descended from each node in nodes
#' @importFrom phangorn allDescendants
#' @keywords internal
#' @export
CladeSizes <- function (tree, nodes) {
  if (is.null(treeOrder <- attr(tree, 'order')) || treeOrder != 'postorder') tree <- Postorder(tree)
  vapply(phangorn::allDescendants(tree)[nodes], length, integer(1))
}
    

#' node depth
#' Wrapper for the ape function
## @useDynLib TreeSearch ape_node_depth
#' @keywords internal
#' @export
C_node_depth <- function (nTip, nNode, parent, child, nEdge) {
  .C("ape_node_depth", as.integer(nTip), as.integer(nNode), as.integer(parent), 
     as.integer(child), as.integer(nEdge), double(nTip + nNode), 1L)[[6]]
}