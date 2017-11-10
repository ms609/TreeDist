#' @title Calculate parsimony score with inapplicable data
#'
#' @description Uses code modified from the Morphy library to calculate a parsimony score 
#' in datasets that contain inapplicable data
#'
#' @template treeParam 
#' @template datasetParam
#' 
#' @examples
#' data(inapplicable.datasets)
#' tree <- TreeSearch::RandomTree(inapplicable.phyData[[1]])
#' result <- InapplicableFitch(tree, inapplicable.phyData[[1]])
#' 
#' @return This function returns the elements from a list containing:
#'    \itemize{
#' \item     The total parsimony score
#' \item     The parsimony score associated with each character 
#' \item     A matrix comprising character reconstructions for each node after the final pass
#'   }
#' The elements to return are specified by the parameter \code{detail}.  
#' If a single element is requested (default) then just that element will be returned
#' If multiple elements are requested then these will be returned in a list.
#' 
#' @seealso \code{\link{MorphyDat}}
#' @seealso \code{\link{BasicSearch}}
#' 
#' @author Martin R. Smith (using C code adapted from MorphyLib, author Martin Brazeau)
#' @importFrom phangorn phyDat
#' @export
InapplicableFitch <- function (tree, dataset) {
  if (class(dataset) != 'phyDat') stop('Invalid data type ', class(dataset), '; should be phyDat.')
  tree <- RenumberTips(Renumber(tree), names(dataset))
  morphyObj <- PhyDat2Morphy(dataset)
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  MorphyTreeLength(tree, morphyObj)
}

#' Calculate parsimony score with inapplicable data
#' 
#' @template labelledTreeParam
#' @template morphyObjParam
#'
#' @return The length of the tree (after weighting)
#'
#' @seealso PhyDat2Morphy
#'
#' @author Martin R. Smith
#' @keywords internal
#' @export
MorphyTreeLength <- function (tree, morphyObj) {
  nTaxa <- mpl_get_numtaxa(morphyObj)
  if (nTaxa != length(tree$tip.label)) stop ("Number of taxa in morphy object (", nTaxa, ") not equal to number of tips in tree")
  treeOrder <- attr(tree, 'order')
  inPostorder <- (!is.null(treeOrder) && treeOrder == "postorder")
  tree.edge <- tree$edge
  # Return:
  ParentChildMorphyLength(tree.edge[, 1], tree.edge[, 2], morphyObj, inPostorder, nTaxa)
}

#' @describeIn MorphyTreeLength Faster function that requires internal tree parameters
#' @template treeParent
#' @template treeChild
#' @author Martin R. Smith
#' @keywords internal
#' @export
ParentChildMorphyLength <- function (parent, child, morphyObj, inPostorder=FALSE, nTaxa=mpl_get_numtaxa(morphyObj)) {
  if (!inPostorder) {
    edgeList <- PostorderEdges(parent, child, nTaxa=nTaxa)
    parent <- edgeList[[1]]
    child <- edgeList[[2]]
  }
  if (nTaxa < 1L) stop("Error: ", mpl_translate_error(nTaxa))
  maxNode <- nTaxa + mpl_get_num_internal_nodes(morphyObj)
  rootNode <- nTaxa + 1L
  allNodes <- rootNode:maxNode
  
  parentOf <- parent[match(1:maxNode, child)]
  # parentOf[rootNode] <- maxNode + 1 # Root node's parent is a dummy node
  parentOf[rootNode] <- rootNode # Root node's parent is a dummy node
  leftChild <- child[length(parent) + 1L - match(allNodes, rev(parent))]
  rightChild <- child[match(allNodes, parent)]
  
  # Return:
  .Call('C_MORPHYLENGTH', as.integer(parentOf -1L), as.integer(leftChild -1L), 
               as.integer(rightChild -1L), morphyObj)
}

#' @describeIn MorphyTreeLength Fastest function that requires internal tree parameters
#' @template parentOfParam
#' @template leftChildParam
#' @template rightChildParam
#' @author Martin R. Smith
#' @keywords internal
#' @export
MorphyLength <- function (parentOf, leftChild, rightChild, morphyObj) {
  # Return:
  .Call('C_MORPHYLENGTH', as.integer(parentOf), as.integer(leftChild), 
               as.integer(rightChild), morphyObj)
}

#' @describeIn MorphyTreeLength Direct call to C function. Use with caution.
#' @param parentOf For each node, numbered in postorder, the number of its parent node.
#' @param leftChild  For each internal node, numbered in postorder, the number of its left 
#'                   child node or tip.
#' @param rightChild For each internal node, numbered in postorder, the number of its right
#'                   child node or tip.
#' @keywords internal
#' @export
C_MorphyLength <- function (parentOf, leftChild, rightChild, morphyObj) {
  .Call('C_MORPHYLENGTH', as.integer(parentOf -1L), as.integer(leftChild -1L), 
               as.integer(rightChild -1L), morphyObj)
}