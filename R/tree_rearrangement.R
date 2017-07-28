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
  as.integer(nEdge))[1:2]), ncol=2)
}

# #' @useDynLib TreeSearch OENN
# RenumberTree <- function (parent, child, nTips, nEdge) {
#   matrix(unlist(.Call('order_edges_number_nodes', as.integer(parent), as.integer(child),
#   as.integer(nTips-1L), as.integer(nEdge))), ncol=2)
# }

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