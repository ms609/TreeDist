
#' Rearrange phylogenetic tree
#' @details \code{RearrangeTree} performs one tree rearrangement of a specified type
#' 
#' @param tree a rooted bifurcating phylogenetic tree with the desired outgroup, with its labels
#'             in an order that matches the Morphy object, and the attributes
#'             \code{pscore}, the tree's parsimony score, and 
#'             \code{hits}, the number of times the best score has been hit in the calling function;
#' @param data a dataset in the format required by ParsimonyScorer
#' @param Rearrange a rearrangement function: probably one of 
#'     \code{\link{RootedNNI}}, \code{\link{RootedSPR}} or \code{\link{RootedTBR}};
#' @param ParsimomyScorer a function that returns a score to be optimised
#' @param  min.score trees longer than \code{min.score}, probably the score of the starting tree,
#'     will be discarded;
#' @template concavityParam 
#' @param  return.single returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE};}
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
#' data('SigSut')
#' random.tree <- RandomTree(SigSut.phy)
#' RearrangeTree(random.tree, SigSut.phy, RootedNNI)
#' 
#' @importFrom parallel clusterCall
#' @export
RearrangeTree <- function (tree, data, Rearrange = NNI, ParsimonyScorer = phangorn::fitch, 
                           min.score=NULL, return.single=TRUE, iter='<unknown>', cluster=NULL,
                           track=0) {
  if (is.null(attr(tree, 'score'))) best.score <- 1e+07 else best.score <- attr(tree, 'score')
  if (is.null(attr(tree, 'hits'))) hits <- 1 else hits <- attr(tree, 'hits')
  if (is.null(cluster)) {
    re.tree <- Rearrange(tree)
    trees <- list(re.tree)
    min.score <- ParsimonyScorer(re.tree, data)
    best.trees <- c(TRUE)
  } else {
    #candidates <- clusterCall(cluster, function(re, tr, k) {ret <- re(tr); attr(ret, 'score') <- ParsimonyScorer(ret, cl.data, k); ret}, Rearrange, tree)
    #scores <- vapply(candidates, function(x) attr(x, 'ps'), 1)
    candidates <- clusterCall(cluster, Rearrange, tree)
    scores <- vapply(candidates, ParsimonyScorer, 1, data, target=min.score) # ~3x faster to do this in serial in r233.
    min.score <- min(scores)
    best.trees <- scores == min.score
    trees <- candidates[best.trees]
  }
  if (best.score < min.score) {
    if (track > 3) cat("\n    . Iteration", iter, '- Min score', min.score, ">", best.score)
  } else if (best.score == min.score) {
    hits <- hits + sum(best.trees)
    if (track > 2) cat("\n    - Iteration", iter, "- Best score", min.score, "hit", hits, "times")
  } else {
    hits <- sum(best.trees)
    if (track > 1) cat("\n    * Iteration", iter, "- New best score", min.score, "found on", hits, "trees")
  }
  if (length(return.single) && return.single) trees <- sample(trees, 1)[[1]]
  attr(trees, 'hits') <- hits
  attr(trees, 'score') <- min.score
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
#' Data(SigSut) # Loads the phyDat object SigSut.phy
#' tree <- RandomTree(SigSut.phy) # 
#' tree <- RenumberTips(tree, names(SigSut.phy))
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
#' @param edge.to.break the index of an edge to bisect, generated randomly if not specified
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
  chosenInternalEdge <- sample(which(child > nTips), 1)
  if(is.na(chosenInternalEdge)) return(NULL)
  nEdge <- length(parent)
  nNode <- tree$Nnode
  if (nNode == 1) return(tree)
  
  end1  <- parent[chosenInternalEdge]
  end2  <- child[chosenInternalEdge]
  
  ind1    <- which(parent == end1)
  ind1    <- ind1[ind1 != chosenInternalEdge][1]
  ind2    <- which(parent == end2)[sample(2, 1)]
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
  ind     <- sample(which(sampleableChild > nTips), 1)
  if(is.na(ind)) return(NULL)
  
  rootNode <- nTips + 1L
  nEdge <- length(parent)
  nNode <- tree$Nnode
  
  if (nNode == 1) return(tree)
  p1      <- parent[ind]
  p2      <- child[ind]
  ind1    <- which(parent == p1)
  ind1    <- ind1[ind1 != ind][1]
  ind2    <- which(parent == p2)[sample(2, 1)]
  new_ind <- c(ind2, ind1)
  old_ind <- c(ind1, ind2)
  child_swap <- child[new_ind]
  edge [old_ind, 2L] <- child_swap
  child[old_ind] <- child_swap

  
  # TODO do all this in a single C function, to minimise number of calls.
  neworder <- NeworderPhylo(nTips, parent, child, nEdge, 2L)
  tree$edge <- OrderEdgesNumberNodes(edge[neworder, 1], edge[neworder, 2], nTips, nEdge)
  tree
  # TODO validate tree length swaps.
}

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
  edge <- tree$edge; parent <- edge[,1L]; child <- edge[,2L]
  root <- nTips + 1L # Assumes fully-resolved bifurcating tree
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
  prune.node <- sample(pruning.candidates, 1)
  moving.subnodes <- c(prune.node, which(DoDescendants(parent, child, nTips, prune.node)))
  moving.nodes <- c(prune.parent <- parent[child==prune.node], moving.subnodes)
  dont.graft.here <- c(moving.nodes, child[parent==prune.parent])
  graft.candidates <- c(root.children[choose.right + 1L], pruning.candidates)
  graft.candidates <- c(graft.candidates[!graft.candidates %in% dont.graft.here])
  graft.node <- sample(graft.candidates, 1)
  graft.edge <- match(graft.node, child)
  graft.parent <- parent[graft.edge]
  graft.child  <-  child[graft.edge]
  
  leading.edge <- match(prune.parent, child)
  prune.edge <- match(prune.node, child)
  parent.duplicate <- parent
  parent.duplicate[prune.edge] <- NA
  sister.edge <- match(prune.parent, parent.duplicate)
  edge[c(leading.edge, sister.edge, graft.edge), 2] <- edge[c(sister.edge, graft.edge, leading.edge), 2]
  
  nEdge <- length(child)
  tree$edge <- OrderEdgesNumberNodes(edge[, 1], edge[, 2], nTips, nEdge)
  tree
}

#' Rooted TBR 
#' @describeIn TBR Perform \acronym{TBR} rearrangement, retaining position of root
#' @importFrom ape is.rooted
#' @importFrom stats runif
#' @export
RootedTBR <- function(tree) { 
  if (!is.rooted(tree)) warning("Tree root is not resolved.  Try:  tree <- SetOutgroup(tree, outgroup).")
  edge <- tree$edge; parent <- edge[,1L]; child <- edge[,2L]
  root <- 1 + (nTips <- dim(edge)[1] - tree$Nnode + 1L)
  root.children <- child[parent==root]
  size <- c(left.size <- sum(DoDescendants(parent, child, nTips, root.children[1L])) + 1L,
            (nTips * 2L - 1L) - left.size - 1L)
  if (min(size) == 1L) {
    outgroup <- tree$tip.label[root.children[size==1L]]
    tree <- TBR(tree)
    return(Root(tree, outgroup))
  } else {
    tip.label <- tree$tip.label
    moves <- (size-3L) * (size-2L) / 2
    subtree.root <- root.children[1L + (runif(1, min=0, max=sum(moves)) > moves[1L])]
    in.crown <- DoDescendants(parent, child, nTips, subtree.root)
    in.crown[subtree.root] <- TRUE
    crown.edges <- parent %in% which(in.crown)
    in.stump <- !in.crown
    in.stump[root] <- FALSE
    stump.edges <- parent %in% which(in.stump)
    stump <- KeepEdges(edge, tip.label, nTips, stump.edges) # faster than DropTip
    crown <- ExtractClade(tree, subtree.root) # faster than KeepEdges
    new.crown <- TBR(crown)
    new.crown$root.edge <- 1
    return (stump + new.crown)
  }
}

#' Subtree Pruning and Rearrangement 
#'
#' Perform one \acronym{SPR} rearrangement on a tree
#' Equivalent to phangorn's kSPR, but 10x faster
#'
#' @template treeParam
#' @author Martin R. Smith
#' @export
SPR <- function(tree) { 
  labels <- tree$tip.label
  nTips  <- length(labels)
  edge   <- tree$edge
  parent <- edge[, 1]
  child  <- edge[, 2]
  nEdge  <- length(child)
  root   <- nTips + 1
  if (nTips < 4) stop ('must be >3 tips for SPR rearrangement!')
  pruning.candidates <- seq(nEdge + 1)[-root]
  repeat {
    prune.node <- sample(pruning.candidates, 1)
    moving.subnodes <- c(prune.node, which(DoDescendants(parent, child, nTips, prune.node)))
    moving.nodes <- c(prune.parent <- parent[child==prune.node], moving.subnodes)
    dont.graft.here <- c(moving.nodes, child[parent==prune.parent])
    graft.node <- c(pruning.candidates[!pruning.candidates %in% dont.graft.here])
    if (length(graft.node) > 1) graft.node <- sample(graft.node, 1)
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

#' TBR
#' 
#' Tree bisection and reconnection
#'
#' \code{TBR} performs a single random \acronym{TBR} iteration.
#'
#' @param tree A bifurcating tree of class \code{\link{phylo}}, with all nodes resolved;
#' @param edge.to.break the index of an edge to bisect, generated randomly if not specified.
#' 
#' @details Branch lengths are not (yet) supported.
#' 
#' @return This function returns a tree in \code{phyDat} format that has undergone one \acronym{TBR} iteration.
#' @references The \acronym{TBR} algorithm is summarized in
#' Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
#' 
#' 
#' @author Martin R. Smith
#' 
#' @seealso RootedTBR useful when the position of the root node should be retained.
#' 
#' @examples{
#' library('ape')
#' tree <- rtree(20, br=NULL)
#' TBR(tree)
#' }
#' @export
TBR <- function(tree, edge.to.break=NULL) { # FROM INAPPLICABLE
# Improvement targets: Root; ExtractClade; DropTip
  nTips <- tree$Nnode + 1
  if (nTips < 3) return (tree)
  tree.edge <- tree$edge
  tree.parent <- tree.edge[, 1]
  tree.child  <- tree.edge[, 2]
  if (nTips == 3) return (Root(tree, sample(tree.child[tree.parent==max(tree.parent)], 1L)))
  all.nodes <- seq_len(2 * (nTips - 1))
  root <- nTips + 1
  if (is.null(edge.to.break)) edge.to.break <- sample(2:nrow(tree.edge), 1) # Only include one root edge
  subtree.root <- tree.child[edge.to.break]
  #cat("\n - ", edge.to.break, subtree.root)
  stump <- if (subtree.root <= nTips) {
    DropTipNoSubtree(tree, subtree.root, root.edge=1)
  } else {
    in.crown <- DoDescendants(tree.parent, tree.child, nTips, subtree.root, just.tips=TRUE)
    DropTipNoSubtree(tree, which(in.crown), root.edge=1)
  }
  stump.len <- dim(stump$edge)[1]
  crown <- ExtractClade(tree, subtree.root) # ~ 2x faster than DropTip
  crown.edge <- crown$edge
  crown.len <- dim(crown.edge)[1]  
  if (crown.len > 1) {
    if (crown.len == 2) {
      rerooted.crown <- crown
    } else {
      crown.parent <- crown.edge[, 1]
      crown.child <- crown.edge[, 2]
      crown.nNode <- crown$Nnode
      crown.tips <- crown.nNode + 1L
      crown.root <- min(crown.parent)
      new.root.candidates <- crown.child[-1] # Include existing root once only
      new.root.node <- sample(new.root.candidates, 1L)
      if (new.root.node <= crown.tips) new.outgroup <- new.root.node else new.outgroup <- which(DoDescendants(crown.parent, crown.child, crown.tips, new.root.node, just.tips=TRUE))
      rerooted.crown <- Root(crown, new.outgroup)
    }
    rerooted.crown$root.edge <- 1L
    if (stump.len > 1) {
      bind.location <- sample(seq_len(stump.len), 1L)
      ret <- BindTree(stump, rerooted.crown, position=1, where=bind.location)
    } else {
      ret <- AddTip(rerooted.crown, crown.root, stump$tip.label)
    }
  } else {
    bind.location <- stump$edge[sample(2L:stump.len, 1L), 2L]
    ret <- AddTip(stump, bind.location, crown$tip.label)
  }
  Renumber(ret)
}
