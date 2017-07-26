#' Rearrange phylogenetic tree
#' @details \code{RearrangeTree} performs one tree rearrangement of a specified type
#' 
#' @param tree a rooted bifurcating phylogenetic tree with the desired outgroup, with its labels
#'             in an order that matches the Morphy object, and the attributes
#'             \code{pscore}, the tree's parsimony score, and 
#'             \code{hits}, the number of times the best score has been hit in the calling function;
#' @param data a dataset in the format required by TreeScorer
#' @param Rearrange a rearrangement function: probably one of 
#'     \code{\link{RootedNNI}}, \code{\link{RootedSPR}} or \code{\link{RootedTBR}};
#' @param TreeScorer a function that returns a score to be optimised
#' @param  minScore trees longer than \code{minScore}, probably the score of the starting tree,
#'     will be discarded;
#' @template concavityParam 
#' @param  returnSingle returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE};}
#'   \item{iter}{iteration number of calling function, for reporting to user only;
#' @param  cluster a cluster, prepared with \code{\link{PrepareCluster}}, to accelerate 
#'     searches on multicore machines;
#' @param verbosity determines how much information to output to screen.
#' 
#' @return{This function returns the most parsimonious of the trees generated, with attributes \code{hits} and \code{pscore}
#'  as described for argument \code{tree}, and with tip labels ordered to match morphyObj.}
#' @author Martin R. Smith
#' @seealso
#'   \itemize{
#'     \item \code{\link{RootedNNI}}
#'     \item \code{\link{RootedSPR}}
#'     \item \code{\link{RootedTBR}}
#'   }
#' 
#' @examples
#' data('Lobo')
#' random.tree <- RandomTree(Lobo.phy)
#' RearrangeTree(random.tree, Lobo.phy, RootedNNI)
#' 
#' @importFrom parallel clusterCall
#' @export
RearrangeTree <- function (tree, data, Rearrange = NNI, TreeScorer = FitchScore, 
                           minScore=NULL, returnSingle=TRUE, iter='<unknown>', cluster=NULL,
                           track=0) {
  if (is.null(attr(tree, 'score'))) bestScore <- 1e+07 else bestScore <- attr(tree, 'score')
  if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
  if (is.null(cluster)) {
    rearrTree <- Rearrange(tree)
    trees <- list(rearrTree)
    minScore <- TreeScorer(rearrTree, data)
    bestTrees <- c(TRUE)
  } else {
    #candidates <- clusterCall(cluster, function(re, tr, k) {ret <- re(tr); attr(ret, 'score') <- TreeScorer(ret, cl.data, k); ret}, Rearrange, tree)
    #scores <- vapply(candidates, function(x) attr(x, 'ps'), 1)
    candidates <- clusterCall(cluster, Rearrange, tree)
    scores <- vapply(candidates, TreeScorer, 1, data, target=minScore) # ~3x faster to do this in serial in r233.
    minScore <- min(scores)
    bestTrees <- scores == minScore
    trees <- candidates[bestTrees]
  }
  if (bestScore < minScore) {
    if (track > 3) cat("\n    . Iteration", iter, '- Min score', minScore, ">", bestScore)
  } else if (bestScore == minScore) {
    hits <- hits + sum(bestTrees)
    if (track > 2) cat("\n    - Iteration", iter, "- Best score", minScore, "hit", hits, "times")
  } else {
    hits <- sum(bestTrees)
    if (track > 1) cat("\n    * Iteration", iter, "- New best score", minScore, "found on", hits, "trees")
  }
  if (length(returnSingle) && returnSingle) trees <- sample(trees, 1L)[[1]]
  attr(trees, 'hits') <- hits
  attr(trees, 'score') <- minScore
  trees
}

#' @useDynLib TreeSearch ape_neworder_phylo
NeworderPhylo <- function (nTaxa, parent, child, nb.edge, whichwise) {
  .C('ape_neworder_phylo', as.integer(nTaxa), as.integer(parent), as.integer(child), 
     as.integer(nb.edge), integer(nb.edge), as.integer(whichwise), NAOK = TRUE)[[5]]
}

#' @useDynLib TreeSearch ape_neworder_pruningwise
NeworderPruningwise <- function (nTaxa, nb.node, parent, child, nb.edge) {
  .C('ape_neworder_pruningwise', as.integer(nTaxa), as.integer(nb.node), as.integer(parent), 
     as.integer(child), as.integer(nb.edge), integer(nb.edge))[[6]]
}

#' @useDynLib TreeSearch order_edges_number_nodes
OrderEdgesNumberNodes <- function (parent, child, nTips, nEdge) {
  matrix(unlist(.C('order_edges_number_nodes', as.integer(parent), as.integer(child),
  as.integer(nTips-1L), as.integer(nEdge))[1:2]), ncol=2)
}

#' Reorder tree Cladewise
#' 
#' A wrapper for \code{ape:::.reorder_ape}
#'
#' @template treeParam
#' @param nTaxa (optional) number of tips in the tree
#' @param edge (optional) the value of tree$edge
#'
#' @return A tree with nodes numbered in postorder
#' @author Modified by Martin R. Smith from \code{.reorder_ape} in \pkg{ape} (Emmanuel Paradis)
#'
#' @keywords internal
#' @export
Cladewise <- function (tree, nTaxa = NULL, edge = tree$edge) {
  if (!is.null(attr(tree, "order")) && attr(tree, "order") == "cladewise") return(tree)
  if (is.null(nTaxa)) nTaxa <- length(tree$tip.label)
  if (is.null(edge)) edge <- tree$edge
  nb.edge <- dim(edge)[1]
  nb.node <- tree$Nnode
  if (nb.node == 1) return(tree)
  if (nb.node >= nTaxa) stop("tree apparently badly conformed")
  
  neworder <- NeworderPhylo(nTaxa, edge[, 1], edge[, 2], nb.edge, 1)
                 
  tree$edge <- edge[neworder, ]
  if (!is.null(tree$edge.length)) tree$edge.length <- tree$edge.length[neworder]
  attr(tree, "order") <- "cladewise"
  tree
}


#' @describeIn Cladewise Reorder tree in Postorder
#' @export
Postorder <- function (tree, nTaxa = length(tree$tip.label), edge = tree$edge) {
  if (!is.null(attr(tree, "order")) && attr(tree, "order") == "postorder") return(tree)
  nb.edge <- dim(edge)[1]
  nb.node <- tree$Nnode
  if (nb.node == 1) return(tree)
  if (nb.node >= nTaxa) stop("tree apparently badly conformed")
  neworder <- NeworderPhylo(nTaxa, edge[, 1], edge[, 2], nb.edge, 2)
  tree$edge <- edge[neworder, ]
  if (!is.null(tree$edge.length)) tree$edge.length <- tree$edge.length[neworder]
  attr(tree, "order") <- "postorder"
  tree
}

#' @describeIn Cladewise Reorder tree Pruningwise
#' @export
Pruningwise <- function (tree, nTaxa = length(tree$tip.label), edge = tree$edge) {
  if (!is.null(attr(tree, "order")) && attr(tree, "order") == 'pruningwise') return(tree)
  nb.edge <- dim(edge)[1]
  nb.node <- tree$Nnode
  if (nb.node == 1) return(tree)
  if (nb.node >= nTaxa) stop("tree apparently badly conformed")
  tree <- Cladewise(tree, nTaxa, edge)
  neworder <- NeworderPruningwise(nTaxa, nb.node, tree$edge[, 1], tree$edge[, 2], nb.edge)
  tree$edge <- tree$edge[neworder, ]
  if (!is.null(tree$edge.length)) tree$edge.length <- tree$edge.length[neworder]
  attr(tree, "order") <- 'pruningwise'
  tree
}

#' Reorder tips
#'
#' \code{RenumberTips(tree, tipOrder)} sorts the tips of a phylogenetic tree 
#' such that the indices in \code{tree$edge[, 2]} correspond to the order of
#' tips given in \code{tipOrder}
#'
#' @template treeParam
#' @param tipOrder A character vector containing the values of 
#'        \code{tree$tip.label} in the desired sort order
#' 
#' @examples
#' Data(Lobo) # Loads the phyDat object Lobo.phy
#' tree <- RandomTree(Lobo.phy) # 
#' tree <- RenumberTips(tree, names(Lobo.phy))
#'
#' @author Martin R. Smith
#' @export
RenumberTips <- function (tree, tipOrder) {
  startOrder <- tree$tip.label
  if (identical(startOrder, tipOrder)) return (tree)
  
  nTip <- length(startOrder)
  child <- tree$edge[, 2]
  tips <- child <= nTip
  
  tree$edge[tips, 2] <- match(startOrder, tipOrder)
  tree$tip.label <- tipOrder
  tree
}


#' @title Tree rearrangement functions
#' 
#' These functions performs a single random \acronym{TBR}, \acronym{SPR} or \acronym{NNI} iteration.
#'
#' Performs a single iteration of the nearest-neigbour interchange, subtree pruning and regrafting,
#' or tree bisection and reconnection algorithms.
#' NNI and SPR are based on the corresponding phangorn functions, but have been re-coded to 
#' improve their speed.
#' 
#' Branch lengths are not supported.
#' 
#' @return Returns a tree with class \code{phylo}.
#'
#' @template treeParam
#' @param edgeToBreak the index of an edge to bisect, generated randomly if not specified
#' 
#' @references
#' The algorithms are summarized in
#' Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
#' 
#' @author Martin R. Smith
#' 
#' @examples
#' tree <- ape:::rtree(20, br=NULL)
#' NNI(tree)
#' SPR(tree)
#' TBR(tree)
#' @export
NNI <- function (tree) {
  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  nTips  <- length(tree$tip.label)
  rootNode <- nTips + 1L
  chosenInternalEdge <- SampleOne(which(child > nTips))
  if(is.na(chosenInternalEdge)) return(NULL)
  nEdge <- length(parent)
  nNode <- tree$Nnode
  if (nNode == 1) return(tree)
  
  end1  <- parent[chosenInternalEdge]
  end2  <- child[chosenInternalEdge]
  ind1    <- which(parent == end1)
  ind1    <- ind1[ind1 != chosenInternalEdge][1]
  ind2    <- which(parent == end2)[sample.int(2L, 1L, useHash=FALSE)]
  new_ind <- c(ind2, ind1)
  old_ind <- c(ind1, ind2)
  child_swap <- child[new_ind]
  child[old_ind] <- child_swap
  tree$edge <- OrderEdgesNumberNodes(parent, child, nTips, nEdge)
  tree
}

#' Rearrange a rooted tree
#'
#' This function performs a rearrangement iteration on a tree, retaining the position of the root.
#'
#' A single \acronym{NNI}, \acronym{SPR} or \acronym{TBR} rearrangement is performed, subject to the constraint that 
#' no taxon may be moved to the opposite side of the root node.
#' Branch lengths are not (yet) supported.
#' 
#' @usage
#' RootedNNI(tree)
#' RootedSPR(tree)
#' RootedTBR(tree)
#'
#' @param tree A bifurcating tree of class \code{\link{phylo}}, with all nodes resolved;
#' 
#' @return This function returns a tree, in \code{phylo} format.
#'
#' @author Martin R. Smith
#' \code{RootedNNI} is abridged from the \pkg{phangorn} function \code{nnin}
#' 
#' @seealso
#' \itemize{
#' \item \code{\link{SetOutgroup}}, set the outgroup of the phylogenetic tree
#' \item \code{\link{NNI}}, unrooted \acronym{NNI} and \acronym{SPR}
#' \item \code{\link{TBR}}, unrooted \acronym{TBR}
#' }
#' 
#' @examples{
#'   require('ape')
#'   tree <- read.tree(text='(((a,b),c),(d,(e,f)));')
#'   tree <- SetOutgroup(tree, c('e', 'f'))
#'   plot(tree)
#'   dev.new()
#'   plot(RootedNNI(tree))
#'   plot(RootedSPR(tree))
#'   plot(RootedTBR(tree))
#' }
#' 
#'
#' @export
RootedNNI <- function (tree) {
  nTips  <- length(tree$tip.label)
  edge    <- tree$edge
  parent  <- edge[, 1]
  child   <- edge[, 2]
  sampleableChild <- child
  sampleableChild[which(parent == as.integer(parent[!match(parent, child, 0)][1]))] <- -1 # Don't want to switch across the root
  chosenInternalEdge <- SampleOne(which(sampleableChild > nTips))
  if(is.na(chosenInternalEdge)) return(NULL)
  
  rootNode <- nTips + 1L
  nEdge <- length(parent)
  nNode <- tree$Nnode
  if (nNode == 1) return(tree)
  
  end1  <- parent[chosenInternalEdge]
  end2  <- child[chosenInternalEdge]
  ind1    <- which(parent == end1)
  ind1    <- ind1[ind1 != chosenInternalEdge][1]
  ind2    <- which(parent == end2)[sample.int(2L, 1L, useHash=FALSE)]
  new_ind <- c(ind2, ind1)
  old_ind <- c(ind1, ind2)
  child_swap <- child[new_ind]
  child[old_ind] <- child_swap
  tree$edge <- OrderEdgesNumberNodes(parent, child, nTips, nEdge)
  tree
}

#' Subtree Pruning and Rearrangement 
#'
#' Perform one \acronym{SPR} rearrangement on a tree
#' 
#' Equivalent to phangorn's kSPR, but faster.
#' Note that rearrangements that only change the position of the root WILL be returned by 
#' \code{SPR}.  If the position of the root is irrelevant (as in Fitch parsimony, for example)
#' then this function will occasionally return a functionally equivalent topology.  
#' \code{RootIrrelevantSPR} will search tree space more efficiently in these cases.
#'
#' @template treeParam
#' @author Martin R. Smith
#' @export
SPR <- function(tree) { 
  tip.label <- tree$tip.label
  nTips  <- length(tip.label)
  if (nTips < 4) stop ('must be >3 tips for SPR rearrangement!')
  edge   <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  nEdge  <- length(child)
  root   <- nTips + 1
  pruning.candidates <- seq(nEdge + 1)[-root]
  repeat {
    prune.node <- SampleOne(pruning.candidates)
    moving.subnodes <- c(prune.node, which(DoDescendants(parent, child, nTips, prune.node)))
    moving.nodes <- c(prune.parent <- parent[child==prune.node], moving.subnodes)
    dont.graft.here <- c(moving.nodes, child[parent==prune.parent])
    graft.node <- c(pruning.candidates[!pruning.candidates %in% dont.graft.here])
    if (length(graft.node) > 1) graft.node <- SampleOne(graft.node)
    if (any(graft.node)) break;
    pruning.candidates <- pruning.candidates[-match(prune.node, pruning.candidates)]
    if (!any(pruning.candidates)) stop('No place to graft pruned tree')
  } 
  
  graft.edge   <- match(graft.node, child)
  graft.parent <- parent[graft.edge]
  graft.child  <-  child[graft.edge]
  prune.edge   <- match(prune.node, child)
  parent.duplicate <- parent
  parent.duplicate[prune.edge] <- NA
  sister.edge  <- match(prune.parent, parent.duplicate)
  if (prune.parent == root) {
    new.root <- child[parent==root]
    new.root <- new.root[new.root != prune.node]
    edge[sister.edge, 2] <- edge[graft.edge, 2]
    edge[graft.edge,  2] <- root
    new.root.spots <- edge==new.root
    edge[edge == root] <- new.root
    edge[new.root.spots] <- root
  } else {
    leading.edge <- match(prune.parent, child)
    edge[c(leading.edge, sister.edge, graft.edge), 2] <- edge[c(sister.edge, graft.edge, leading.edge), 2]
  }
  
  tree$edge <- OrderEdgesNumberNodes(edge[, 1], edge[, 2], nTips, nEdge)
  tree
}

#' @describeIn SPR Root position irrelevant SPR
#' @export
RootIrrelevantSPR <- function(tree) { 
  tip.label <- tree$tip.label
  nTips  <- length(tip.label)
  if (nTips < 4) stop ('must be >3 tips for SPR rearrangement!')
  edge   <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  nEdge  <- length(child)
  root   <- nTips + 1
  pruning.candidates <- seq(nEdge + 1)[-root]
  repeat {
    prune.node <- SampleOne(pruning.candidates)
    moving.subnodes <- c(prune.node, which(DoDescendants(parent, child, nTips, prune.node)))
    moving.nodes <- c(prune.parent <- parent[child==prune.node], moving.subnodes)
    if (prune.parent == root) {
      rootChildren <- child[parent==root]
      dont.graft.here <- c(moving.nodes, rootChildren, child[parent %in% rootChildren])
    } else {
      dont.graft.here <- c(moving.nodes, child[parent==prune.parent])
    }
    graft.node <- c(pruning.candidates[!pruning.candidates %in% dont.graft.here])
    if (length(graft.node) > 1) graft.node <- SampleOne(graft.node)
    if (any(graft.node)) break;
    pruning.candidates <- pruning.candidates[-match(prune.node, pruning.candidates)]
    if (!any(pruning.candidates)) stop('No place to graft pruned tree')
  } 
  
  graft.edge   <- match(graft.node, child)
  graft.parent <- parent[graft.edge]
  graft.child  <-  child[graft.edge]
  prune.edge   <- match(prune.node, child)
  parent.duplicate <- parent
  parent.duplicate[prune.edge] <- NA
  sister.edge  <- match(prune.parent, parent.duplicate)
  if (prune.parent == root) {
    new.root <- child[parent==root]
    new.root <- new.root[new.root != prune.node]
    edge[sister.edge, 2] <- edge[graft.edge, 2]
    edge[graft.edge,  2] <- root
    new.root.spots <- edge==new.root
    edge[edge == root] <- new.root
    edge[new.root.spots] <- root
  } else {
    leading.edge <- match(prune.parent, child)
    edge[c(leading.edge, sister.edge, graft.edge), 2] <- edge[c(sister.edge, graft.edge, leading.edge), 2]
  }
  
  tree$edge <- OrderEdgesNumberNodes(edge[, 1], edge[, 2], nTips, nEdge)
  tree
}

# interestingTree <- read.tree(text='(((a, (b, c)), (d, e)), (f, g));')
# interestingPruneNode <- 9
# interestingGraftNode <- 6

#' Rooted SPR rearrangement
#'
#' @importFrom ape is.rooted 
#' @importFrom stats runif 
#' @describeIn SPR Perform \acronym{SPR} operation, retaining position of root
#' @export
RootedSPR <- function(tree) {
  if (!is.rooted(tree)) warning("Tree root is not resolved.  Try:  tree <- SetOutgroup(tree, outgroup).")
  tip.label <- tree$tip.label
  nTips <- length(tip.label)
  if (nTips < 4) stop("Must be >3 tips for SPR rearrangement.")
  edge <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  nEdge  <- length(child)
  root   <- nTips + 1
  
  root.children <- child[parent==root]
  left.nodes <- DoDescendants(parent, child, nTips, root.children[1L])
  right.nodes <- !left.nodes
  left.nodes[root.children] <- right.nodes[root.children] <- right.nodes[root] <- FALSE
  size <- c(sum(left.nodes), sum(right.nodes))
  moves <- (size-2L) * (size-1L) / 2
  moves[size < 3] <- 0
  if (!max(moves)) return (tree)
  
  choose.right <- runif(1, min=0, max=sum(moves)) > moves[1]
  pruning.candidates <- if (choose.right) which(right.nodes) else which(left.nodes)
  subtree.base <- child[parent==root.children[choose.right + 1L]]
  subtree.basal.tip <- subtree.base < root
  if (any(subtree.basal.tip)) pruning.candidates <- pruning.candidates[-match(subtree.base[!subtree.basal.tip], pruning.candidates)]
  
  prune.node <- SampleOne(pruning.candidates)
  moving.subnodes <- c(prune.node, which(DoDescendants(parent, child, nTips, prune.node)))
  moving.nodes <- c(prune.parent <- parent[child==prune.node], moving.subnodes)
  dont.graft.here <- c(moving.nodes, child[parent==prune.parent])
  graft.candidates <- c(root.children[choose.right + 1L], pruning.candidates)
  graft.candidates <- graft.candidates[!graft.candidates %in% dont.graft.here]

  graft.child  <- SampleOne(graft.candidates)
  graft.edge   <- match(graft.child, child)
  graft.parent <- parent[graft.edge]
  
  leading.edge <- match(prune.parent, child)
  prune.edge <- match(prune.node, child)
  parent.duplicate <- parent
  parent.duplicate[prune.edge] <- NA
  sister.edge <- match(prune.parent, parent.duplicate)
  edge[c(leading.edge, sister.edge, graft.edge), 2] <- edge[c(sister.edge, graft.edge, leading.edge), 2]
  
  tree$edge <- OrderEdgesNumberNodes(edge[, 1], edge[, 2], nTips, nEdge)
  tree
}

