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
#'  Single taxon tree
#'
#' Create a phylogenetic 'tree' that comprises a single taxon.
#'
#' @usage SingleTaxonTree(label)
#' @param   label a character vector specifying the label of the tip.
#' @return This function returns a \code{phylo} object containing a single tip with the specified label.
#' @seealso \code{\link{TwoTipTree}}
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
#' library('phangorn')
#' tree <- rtree(20, br=NULL)
#' plot(tree); nodelabels(); nodelabels(33, 33, bg='yellow'); dev.new()
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

#' SetOutgroup
#' 
#'  Root a phylogenetic tree
#'
#' Sets the root of a phylogenetic tree, such that one child of the root node is \code{outgroup}.
#'
#' @usage SetOutgroup(tree, outgroup)
#' 
#' @param tree A tree of class \code{\link{phylo}}, with all nodes resolved;
#' @param outgroup a vector of mode numeric or character specifying the new outgroup.
#' 
#' @return This function returns a rooted tree with the new outgroup.
#' 
#' @author Martin R. Smith
#' 
#' 
#' @seealso root ought to produce the same result as this function,
#'  but does not always do so in practice.
#' 
#' @examples{
#'   require('ape')
#'   tree <- read.tree(text='(((a,b),c),(d,e));')
#'   plot(tree)
#'   plot(Root(tree, c('a', 'b')))
#'   plot(Root(tree, 3))
#' }
#' 
#' @aliases SetOutgroup
#' @aliases Root
#' 
#' @export SetOutgroup Root
SetOutgroup <- function (tree, outgroup) {
  if (class(tree) != 'phylo') stop ('"tree" must be of class "phylo"')
  tip <- tree$tip.label
  if (is.character(outgroup)) {outgroup <- match(outgroup, tip, nomatch=0); outgroup <- outgroup[as.logical(outgroup)]}
  if (length(outgroup) < 1) stop ('"outgroup" not specified')
  if (!is.null(tree$edge.length)) {tree$edge.length <- NULL; warning('Edge lengths are not supported and have been dropped.')}
  nTips <- length(tip)
  edge <- tree$edge
  parent <- edge[,1]
  child <- edge[,2]
  root <- min(parent)
  root.children <- child[parent==root]
  rooted.on.ingroup <- FALSE
  # Check that outgroup is currently monophyletic, swapping with 'ingroup' if it's not
  repeat { 
    ancestry <- GetAncestors(parent, child, outgroup)
    if (length(outgroup) > 1) {
      common.ancestors <- Reduce(intersect, ancestry)
      outgroup.root.node <- max(common.ancestors)
    } else {
      common.ancestors <- c(outgroup, ancestry)
      outgroup.root.node <- outgroup
    }
    if (outgroup.root.node != root) break
    if (rooted.on.ingroup) stop ('Cannot root tree: polyphyletic outgroup straddles root')
    outgroup <- seq_along(tip)[-outgroup] # outgroup straddles root; root on ingroup instead
    rooted.on.ingroup <- TRUE
  }
 
  visit.node.forwards <- function (old.tree.node, parent.number, new.edges) {
    blank.edge <- which.min(new.edges) # number of first as-yet-unspecified edge
    if (old.tree.node <= nTips) { # Adding a tip
      new.edges[blank.edge, ] <- c(parent.number, old.tree.node)
    } else { # Adding a node
      this.node.new.number <- max(c(nTips + 1, new.edges[, 2])) + 1
      new.edges[blank.edge, ] <- c(parent.number, this.node.new.number)    
      descendant.nodes <- DoDescendants(parent, child, nTips, old.tree.node)
      n.new.nodes <- sum(descendant.nodes)
      fill.edge.index <- blank.edge + (1:n.new.nodes)
      fill.edge <- edge[child %in% which(descendant.nodes), ]
      node.spots <- fill.edge > nTips
      fill.edge[node.spots] <- fill.edge[node.spots] + this.node.new.number - old.tree.node
      new.edges[fill.edge.index, ] <- fill.edge
    }
    new.edges
  }
  visit.node.backwards <- function (arrival.edge, last.node.number, new.edges) {
    previous.node <- child[arrival.edge]
    this.node <- parent[arrival.edge]
    forward.node <- child[parent==this.node]
    forward.node <- forward.node[forward.node != previous.node]
    blank.edge <- which.min(new.edges) # First edge left blank after visit.node.forwards
    this.node.new.number <- max(c(nTips + 1L, new.edges[,2L])) + 1L
    if (this.node != root) {
      new.edges[blank.edge, ] <- c(last.node.number, this.node.new.number)
      if (length(forward.node) > 0) {
        for (fwd in forward.node) {
          new.edges <- visit.node.forwards(fwd, this.node.new.number, new.edges)
        }
      }
    } else { # Root edge; don't create a new node
      if (length(forward.node) > 0) {
        for (fwd in forward.node) {
          new.edges <- visit.node.forwards(fwd, last.node.number, new.edges)
        }
      }
    }
    backward.edge <- child.index[this.node]
#   cat("\n, this.node", this.node, "backward.edge", backward.edge, "match", match(this.node, child))  # For debugging only
#   arrival.edge <- backward.edge; last.node.number <- this.node.new.number              # For debugging only
    if (!is.na(backward.edge)) new.edges <- visit.node.backwards(backward.edge, this.node.new.number, new.edges)
    new.edges
  }
  
  last.edge <- length(parent)
  child.index <- order(c(child, root))
  child.index[root] <- NA
  new.edges <- matrix(0, last.edge, 2)
  new.edges <- visit.node.forwards(outgroup.root.node, root, new.edges)
  arrival.edge <- child.index[outgroup.root.node]; last.node.number <- root # For debugging only - DELETE
  new.edges <- visit.node.backwards(child.index[outgroup.root.node], root, new.edges)
  tree$edge <- new.edges
  tree
}
#' @export
Root <- SetOutgroup

#' Get Ancestors
#'
#' \code{GetAncestors} gets the ancestors of each node in a tree
#' It's a more efficient version of \code{\link[phangorn]{Ancestors}}
#'
#' @param PARAM is a parameter you should send to it
#' 
#' @examples
#' to_do <- TRUE
#' 
#' @return This function returns :
#'   
#' @author Martin R. Smith
#' @export
GetAncestors <- function (parent, child, node) {
  if (length(node) == 1) {
    pvector <- numeric(max(parent))
    pvector[child] <- parent
    anc <- function(pvector, node) {
      res <- numeric(0)
      repeat {
        anc <- pvector[node]
        if (anc == 0) 
            break
        res <- c(res, anc)
        node <- anc
      }
      res
    }
    return(anc(pvector, node))
  } else AllAncestors(parent, child)[node]
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

#' Get Descendants
#'
#' \code{GetDescendants} returns the descendants of a specified node
#'
#' @template treeParam
#' @param node the number of the internal node whose descendants should be returned
#' @param just.tips , should return value include all nodes or just tips?
#' 
#' @examples
#' tr <- ape::rtree(20)
#' GetDescendants(tr, 25)
#' 
#' @return This function returnsa vector containing descendant nodes in numerical order
#'   
#' @author Martin R. Smith
#' @export
GetDescendants <- function (tree, node, ...) {
  nTip <- length(tree$tip.label)
  edge <- tree$edge
  edge1 <- edge[,1]
  edge2 <- edge[,2]
  return (which(DoDescendants(edge1, edge2, nTip, node, ...)))
}

#' TITLE GOES HERE
#'
#' \code{DoDescendants} does something useful
#'
#' @param edge1 parent nodes: from tree$edge[,1]
#' @param edge2 parent nodes: from tree$edge[,2]
#' @param node  number of an internal node
#' @param just.tips (logical) should return value include all nodes or just tips?
#' 
#' @examples
#' warning(to_do <- TRUE)
#' 
#' @return This function returns a vector containing descendant nodes in numerical order
#'   
#' @author Martin R. Smith
#' @export
DoDescendants <- function (edge1, edge2, nTip = length(edge1) / 2 + 1, node, 
                           just.tips = FALSE, just.internal=FALSE, 
                           include.ancestor = FALSE) {
  is.descendant <- logical((nTip * 2) - 1)
  if (include.ancestor) is.descendant[node] <- TRUE;
  node.children <- function (node, is.descendant) {
    nc <- edge2[edge1 %in% node]
    is.descendant[nc] <- TRUE
    if (length(nc)) is.descendant <- node.children(nc, is.descendant)
    is.descendant
  }
  is.descendant <- node.children(node, is.descendant)
  if (just.tips) return (is.descendant[seq_len(nTip)]) else if (just.internal) is.descendant[seq_len(nTip)] <- FALSE 
  return (is.descendant)
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
  vapply(phangorn:::allDescendants(tree)[nodes], length, integer(1))
}
    

#' TwoTipTree
#'
#' Two-tipped tree
#'
#' This function generates a tree of class \code{\link{phylo}} containing two tips.
#'
#' @usage TwoTipTree(tip1, tip2)
#' 
#' @param tip1 A character string representing the label for tip 1
#' @param tip2 A character string representing the label for tip 2
#' 
#' @return This function returns a \code{phylo} object with a single root node and two tips with
#' the specified labels.
#'
#' @author Martin R. Smith
#' @seealso SingleTaxonTree
#' @seealso bind.tree
#' @examples TwoTipTree('Homo', 'Pan')
#' @keywords  tree 
#' 
#' @export
TwoTipTree <- function (tip1, tip2) {
  ret <- list(
    edge = matrix(c(3L, 3L, 1L, 2L), 2, 2),
    tip.label = c(tip1, tip2),
    Nnode = 1L
  )
  attr(ret, 'class') <- 'phylo'
  attr(ret, 'order') <- 'preorder'
  ret
}

#' Bind trees
#'
#' Copied from \code{ape:::bind.tree}.
#' Changes: 
#'   - use (x|y).edge in place of (x|y)$edge for efficiency
#'   - Delete interactive option
#'   - Use Cladewise instead of reorder (..., 'cladewise')
#'
#' @param x First tree fragment
#' @param y Second tree fragment 
#' @param where which edge to use to bind fragments [?]
#' @param position Where on second tree to bind [?]
#' 
#' @return A single tree comprising the two fragments
#'   
#' @author Martin R. Smith, based on APE (Emmanuel Paradis)
#' @export
BindTree <- function(x, y, where = "root", position = 0) {
    nx <- length(x$tip.label)
    mx <- x$Nnode
    ROOTx <- nx + 1L
    ny <- length(y$tip.label)
    my <- y$Nnode

    
    if (where == 0 || where == "root") where <- ROOTx
    if (position < 0) position <- 0
    if (where > nx + mx)
        stop("argument 'where' out of range for tree 'x'")
    

    ## check whether both trees have branch lengths:
    switch(is.null(x$edge.length) + is.null(y$edge.length) + 1L,
           wbl <- TRUE, {
               x$edge.length <- y$edge.length <- NULL
               wbl <- FALSE
               warning("one tree has no branch lengths, they have been ignored")
           },
           wbl <- FALSE)

    yHasNoRootEdge <- is.null(y$root.edge)
    xHasNoRootEdge <- is.null(x$root.edge)

    x.edge <- x$edge
    y.edge <- y$edge
    ## find the row of 'where' before renumbering
    if (where == ROOTx) case <- 1 else {
        i <- which(x.edge[, 2] == where)
        case <- if (where <= nx) 2 else 3
    }
    ## case = 1 -> y is bound on the root of x
    ## case = 2 -> y is bound on a tip of x
    ## case = 3 -> y is bound on a node of x

    ## check that 'position' is correct
    if (position && wbl) {
### New in ape 3.0-1: this makes possible binding 'y' below
### a node of 'x' thus creating a new node in 'x'
###        if (!wbl)
###            stop("'position' is non-null but trees have no branch lengths")
        if (case == 1) {
            if (xHasNoRootEdge)
                stop("tree 'x' has no root edge")
            if (position > x$root.edge)
                stop("'position' is larger than x's root edge")
        } else {
            if (x$edge.length[i] < position)
                stop("'position' is larger than the branch length")
        }
    }

    ## the special case of substituting two tips:
    if (case == 2 && ny == 1 && !position) {
        x$tip.label[x.edge[i, 2]] <- y$tip.label
        if (wbl)
            x$edge.length[i] <- x$edge.length[i] + y$edge.length
        return(x)
    }

    x <- Cladewise(x)
    y <- Cladewise(y)

### because in all situations internal nodes need to be
### renumbered, they are changed to negatives first, and
### nodes eventually added will be numbered sequentially

    nodes <- x.edge > nx
    x.edge[nodes] <- -(x.edge[nodes] - nx) # -1, ..., -mx
    nodes <- y.edge > ny
    y.edge[nodes] <- -(y.edge[nodes] - ny + mx) # -(mx+1), ..., -(mx+my)
    ROOT <- -1L # may change later
    next.node <- -(mx + my) - 1L

    ## renumber now the tips in y:
    new.nx <- if (where <= nx && !position) nx - 1L else nx
    y.edge[!nodes] <- y.edge[!nodes] + new.nx

    ## if 'y' as a root edge, use it:
    if (!yHasNoRootEdge) {
        y.edge <- rbind(c(0, y.edge[1]), y.edge)
        ##                ^ will be filled later
        next.node <- next.node - 1L
        if (wbl) y$edge.length <- c(y$root.edge, y$edge.length)
    }

    switch(case, { # case = 1
        if (position) {
            x$root.edge <- x$root.edge - position
            x.edge <- rbind(c(next.node, x.edge[1]), x.edge)
            ROOT <- next.node
            if (wbl) x$edge.length <- c(position, x$edge.length)
        }
        if (yHasNoRootEdge) {
            j <- which(y.edge[, 1] == y.edge[1])
            y.edge[j, 1] <- ROOT
        } else y.edge[1] <- ROOT
        x.edge <- rbind(x.edge, y.edge)
        if (wbl)
            x$edge.length <- c(x$edge.length, y$edge.length)
    }, { # case = 2
        if (position) {
            x.edge[i, 2] <- next.node
            x.edge <- rbind(x.edge[1:i, ], c(next.node, where), x.edge[-(1:i), ])
            if (wbl) {
                x$edge.length[i] <- x$edge.length[i] - position
                x$edge.length <- c(x$edge.length[1:i], position, x$edge.length[-(1:i)])
            }
            i <- i + 1L
            if (yHasNoRootEdge) {
                j <- which(y.edge[, 1] == y.edge[1])
                y.edge[j, 1] <- x.edge[i, 1]
            } else y.edge[1] <- x.edge[i, 1]
        } else {
            if (yHasNoRootEdge) x.edge[i, 2] <- y.edge[1]
            else {
                ## the root edge of y is fused with the terminal edge of x
                if (wbl) y$edge.length[1] <- y$edge.length[1] + x$edge.length[i]
                y.edge[1] <- x.edge[i, 1]
                ## delete i-th edge in x:
                x.edge <- x.edge[-i, ]
                if (wbl) x$edge.length <- x$edge.length[-i]
                i <- i - 1L
            }
            x$tip.label <- x$tip.label[-where]
            ## renumber the tips that need to:
            ii <- which(x.edge[, 2] > where & x.edge[, 2] <= nx)
            x.edge[ii, 2] <- x.edge[ii, 2] - 1L
        }
        x.edge <- rbind(x.edge[1:i, ], y.edge, x.edge[-(1:i), ])
        if (wbl)
            x$edge.length <- c(x$edge.length[1:i], y$edge.length, x$edge.length[-(1:i)])
    }, { # case = 3
        if (position) {
            if (yHasNoRootEdge) {
                j <- which(y.edge[, 1] == y.edge[1])
                y.edge[j, 1] <- next.node
            } else y.edge[1] <- next.node
            x.edge <- rbind(x.edge[1:i, ], c(next.node, x.edge[i, 2]), x.edge[-(1:i), ])
            x.edge[i, 2] <- next.node
            if (wbl) {
                x$edge.length[i] <- x$edge.length[i] - position
                x$edge.length <- c(x$edge.length[1:i], position, x$edge.length[-(1:i)])
            }
            i <- i + 1L
        } else {
            if (yHasNoRootEdge) {
                j <- which(y.edge[, 1] == y.edge[1])
                y.edge[j, 1] <- x.edge[i, 2]
            } else y.edge[1] <- x.edge[i, 2]
        }
        x.edge <- rbind(x.edge[1:i, ], y.edge, x.edge[-(1:i), ])
        if (wbl)
            x$edge.length <- c(x$edge.length[1:i], y$edge.length, x$edge.length[-(1:i)])
    })

    x$tip.label <- c(x$tip.label, y$tip.label)

    if (is.null(x$node.label)) {
        if (!is.null(y$node.label))
            x$node.label <- c(rep(NA, mx), y$node.label)
    } else {
        x$node.label <-
            if (is.null(y$node.label)) c(x$node.label, rep(NA, my))
            else c(x$node.label, y$node.label)
    }

    n <- length(x$tip.label)
    x$Nnode <- dim(x.edge)[1] + 1L - n

    ## update the node labels before renumbering (this adds NA for
    ## the added nodes, and drops the label for those deleted)
    if (!is.null(x$node.label))
        x$node.label <- x$node.label[sort(-unique(x.edge[, 1]))]

    ## renumber nodes:
    newNb <- integer(x$Nnode)
    newNb[-ROOT] <- n + 1L
    sndcol <- x.edge[, 2] < 0
    ## executed from right to left, so newNb is modified before x.edge:
    x.edge[sndcol, 2] <- newNb[-x.edge[sndcol, 2]] <- n + (2:x$Nnode)
    x.edge[, 1] <- newNb[-x.edge[, 1]]
    x$edge <- x.edge
    if (!is.null(x$node.label))
        x$node.label <- x$node.label[order(newNb[newNb > 0])]

    x
}

#' node depth
#' Wrapper for the ape function
## @useDynLib TreeSearch ape_node_depth
#' @keywords internal
#' @export
C_node_depth <- function (nTip, nNode, parent, child, nEdge) {
  .C("node_depth", as.integer(nTip), as.integer(nNode), as.integer(parent), 
     as.integer(child), as.integer(nEdge), double(nTip + nNode), 1L, PACKAGE='ape')[[6]]
}