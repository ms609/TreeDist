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

#' Neighbour Joining Tree
#' 
#' Generates a rooted neighbour joining tree, with no edge lengths
#'
#' @template datasetParam
#' 
#' @return an object of class \code{phylo}
#'
#' @author Martin R. Smith
#' @importFrom ape nj root
#' @importFrom phangorn dist.hamming
#' @export
NJTree <- function (dataset) {
  nj.tree <- nj(dist.hamming(dataset))
  nj.tree <- root(nj.tree, outgroup=names(dataset)[1], resolve.root=TRUE)
  nj.tree$edge.length <- NULL
  nj.tree
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
