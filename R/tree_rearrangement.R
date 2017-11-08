#' Rearrange phylogenetic tree
#' @details \code{RearrangeTree} performs one tree rearrangement of a specified type
#' 
#' @param tree a rooted bifurcating phylogenetic tree with the desired outgroup, with its labels
#'             in an order that matches the Morphy object, and the attributes
#'             \code{pscore}, the tree's parsimony score, and 
#'             \code{hits}, the number of times the best score has been hit in the calling function;
#' @template treeScorerParam
#' @param Rearrange a rearrangement function: probably one of 
#'     \code{\link{RootedNNI}}, \code{\link{RootedSPR}} or \code{\link{RootedTBR}};
#' @param  minScore trees longer than \code{minScore}, probably the score of the starting tree,
#'     will be discarded;
#' @param  returnSingle returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE};}
#'   \item{iter}{iteration number of calling function, for reporting to user only;
#' @template verbosityParam
#'
#' @template treeScorerDots
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
#' RearrangeTree(random.tree, dataset=Lobo.phy, Rearrange=RootedNNI)
#' 
#' @export
RearrangeTree <- function (tree, TreeScorer = FitchScore, Rearrange = RootedNNI, 
                           minScore=NULL, returnSingle=TRUE, iter='<unknown>', verbosity=0, ...) {
  if (is.null(attr(tree, 'score'))) bestScore <- 1e+07 else bestScore <- attr(tree, 'score')
  if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
  rearrTree <- Rearrange(tree)
  trees <- list(rearrTree)
  minScore <- TreeScorer(tree=rearrTree, ...)
  bestTrees <- c(TRUE)

  if (bestScore < minScore) {
    if (verbosity > 3) cat("\n    . Iteration", iter, '- Min score', minScore, ">", bestScore)
  } else if (bestScore == minScore) {
    hits <- hits + sum(bestTrees)
    if (verbosity > 2) cat("\n    - Iteration", iter, "- Best score", minScore, "hit", hits, "times")
  } else {
    hits <- sum(bestTrees)
    if (verbosity > 1) cat("\n    * Iteration", iter, "- New best score", minScore, "found on", hits, "trees")
  }
  if (length(returnSingle) && returnSingle) {
    trees <- if (length(trees) > 1) sample(trees, 1L)[[1]] else trees[[1]]
  }
  attr(trees, 'hits') <- hits
  attr(trees, 'score') <- minScore
  trees
}

#' neworder_phylo
#' Wrapper for the ape function
## @useDynLib TreeSearch ape_neworder_phylo
#' @keywords internal
#' @export
NeworderPhylo <- function (nTaxa, parent, child, nb.edge, whichwise) {
  .C(C_ape_neworder_phylo, as.integer(nTaxa), as.integer(parent), as.integer(child), 
     as.integer(nb.edge), integer(nb.edge), as.integer(whichwise), NAOK = TRUE)[[5]]
}

#' neworder_pruningwise
#' Wrapper for the ape function
## @useDynLib TreeSearch ape_neworder_pruningwise
#' @keywords internal
#' @export
NeworderPruningwise <- function (nTaxa, nb.node, parent, child, nb.edge) {
  .C(C_ape_neworder_pruningwise, as.integer(nTaxa), as.integer(nb.node), as.integer(parent), 
     as.integer(child), as.integer(nb.edge), integer(nb.edge))[[6]]
}


#' Order edges and number nodes
#' Wrapper for the C function
#' @return an edge matrix for a tree following the usual convention for edge and node numbering
## @useDynLib TreeSearch order_edges_number_nodes
#' @keywords internal
#' @export
OrderEdgesNumberNodes <- function (parent, child, nTips, nEdge = length(parent)) {
  matrix(unlist(.C(C_order_edges_number_nodes, as.integer(parent), as.integer(child),
  as.integer(nEdge))[1:2]), ncol=2)
}

#' Renumber tree
#' Order edges and number nodes
#' Wrapper for the C function RENUMBER_TREE
#' @return an edge matrix for a tree in following the usual preorder convention for edge and node numbering 
## @useDynLib TreeSearch RENUMBER_TREE
#' @keywords internal
#' @export
RenumberTree <- function (parent, child, nEdge = length(parent)) {
  matrix(.Call(C_RENUMBER_TREE, as.integer(parent), as.integer(child), as.integer(nEdge)), ncol=2)
}

#' @describeIn RenumberTree Instead returns a list containing two items corresponding to the new parent and child vectors
## @useDynLib TreeSearch RENUMBER_EDGES
#' @keywords internal
#' @export
RenumberEdges <- function (parent, child, nEdge = length(parent)) {
  .Call(C_RENUMBER_EDGES, as.integer(parent), as.integer(child), as.integer(nEdge))
}

#' Reorder tree Cladewise
#' 
#' A wrapper for \code{ape:::.reorder_ape}.  Calling this C function directly is approximately twice as fast as using
#' \code{ape::\link[ape]{cladewise}} or \code{ape::\link[ape]{postorder}}
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

#' @describeIn Cladewise Reorder parent and child edges in Postorder
#' @template treeParent
#' @template treeChild
#' @template nTaxaParam
#' @export
PostorderEdges <- function (parent, child, 
                            nEdge = length(parent), 
                            nNode = nEdge / 2L, 
                            nTaxa = nNode + 1L) {
  if (nNode == 1) return(list(parent, child))
  newOrder <- NeworderPhylo(nTaxa, parent, child, nEdge, 2)
  list(parent[newOrder], child[newOrder])
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

#' @describeIn Cladewise Reorder tree in Preorder (special case of cladewise)
#' @export
Preorder <- function (tree) {
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  tree$edge <- RenumberTree(parent, child)
  attr(tree, 'order') <- 'preorder'
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
#' data(Lobo) # Loads the phyDat object Lobo.phy
#' tree <- RandomTree(Lobo.phy)
#' tree <- RenumberTips(tree, names(Lobo.phy))
#'
#' @author Martin R. Smith
#' @export
RenumberTips <- function (tree, tipOrder) {
  startOrder <- tree$tip.label
  if (identical(startOrder, tipOrder)) return (tree)
  if (length(startOrder) != length(tipOrder)) stop("Tree labels and tipOrder must match")
  
  nTip <- length(startOrder)
  child <- tree$edge[, 2]
  tips <- child <= nTip
  
  matchOrder <- match(startOrder, tipOrder)
  if (any(is.na(matchOrder))) stop("All tree labels must occur in tipOrder")  
  tree$edge[tips, 2] <- matchOrder[tree$edge[tips, 2]]
  tree$tip.label <- tipOrder
  tree
}