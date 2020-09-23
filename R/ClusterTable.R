#' Convert phylogenetic tree to `ClusterTable`
#'
#' `as.ClusterTable()` converts a phylogenetic tree to a `ClusterTable` object,
#' which is an internal representation of its splits suitable for rapid tree
#' distance calculation (per Day, 1985).
#'
#' @param x Object to convert into `ClusterTable`: perhaps a tree of class
#'  \code{\link[ape:read.tree]{phylo}}.
#' @param tipLabels Character vector specifying sequence in which to order
#' tip labels.
#' @param \dots Presently unused.
#'
#' @return `as.ClusterTable()` returns an object of class `ClusterTable`.
#'
#' @references \insertRef{Day1985}{TreeDist}
#' @examples
#' tree1 <- ape::read.tree(text = "(A, (B, (C, (D, E))));");
#' tree2 <- ape::read.tree(text = "(A, (B, (D, (C, E))));");
#' ct1 <- as.ClusterTable(tree1)
#' summary(ct1)
#' as.matrix(ct1)
#' 
#' # Tip label order must match ct1 to allow comparison
#' ct2 <- as.ClusterTable(tree2, tipLabels = LETTERS[1:5])
#' 
#' #TODO: RF (ct1, ct2)
#' @template MRS
#' @name ClusterTable
NULL

#' @rdname ClusterTable
#' @export
as.ClusterTable <- function (x, tipLabels = NULL, ...) UseMethod('as.ClusterTable')

#' @rdname ClusterTable
#' @importFrom TreeTools NTip RenumberTips
#' @export
as.ClusterTable.phylo <- function (x, tipLabels = NULL, ...) {
  x <- Preorder(x)
  if (is.null(tipLabels)) {
    tipLabels <- x$tip.label
  } else {
    x <- RenumberTips(x, tipLabels)
  }
  structure(ClusterTable_new(x),
            nTip = NTip(x),
            tip.label = tipLabels,
            class = 'ClusterTable')
}

#' @rdname ClusterTable
#' @export
as.matrix.ClusterTable <- function (x, ...) {
  ClusterTable_matrix(x)
}

#' Print `ClusterTable` object
#'
#' S3 method for objects of class `ClusterTable`.
#'
#' @param x Object of class `ClusterTable`.
#' @param \dots Additional arguments for consistency with S3 method (unused).
#'
#' @export
print.ClusterTable <- function (x, ...) {
  nTip <- attr(x, 'nTip')
  labels <- attr(x, 'tip.label')
  cat("ClusterTable on" , nTip, "leaves:", labels[1], "..", labels[nTip])
}

#' @inherit summary
#' @export
summary.ClusterTable <- function (x, ...) {
  nTip <- attr(x, 'nTip')
  mat <- ClusterTable_matrix(x)
  cat("ClusterTable on" , nTip, "leaves:\n")
  cat(" ", rep(c(1:9, ' '), length.out = nTip), "\n", sep = '')
  apply(mat, 1, function (x) {
    if (x[1] > 0) {
      cat(' ', rep('.', x[1] - 1), rep('*', 1 + x[2] - x[1]), rep('.', nTip - x[2]), "\n", sep = '')
    }
  })
  
  cat(paste0(" ", seq_len(nTip), ": ", attr(x, 'tip.label')[ClusterTable_decode(x)]), "\n")
  
}

#' @importFrom TreeTools TipLabels
#' @inherit TipLabels
#' @export
TipLabels.ClusterTable <- function (x, single = TRUE) attr(x, 'tip.label')

#' @inherit NTip
#' @importFrom TreeTools NTip
#' @export
NTip.ClusterTable <- function (phy) attr(phy, 'nTip')

#' @inherit NSplits
#' @importFrom TreeTools NSplits
#' @export
NSplits.ClusterTable <- function (x) nrow(as.matrix(x)) - 3L # Root + Ingroup + All-leaves

#' @rdname SplitsInBinaryTree
#' @inherit SplitsInBinaryTree
#' @importFrom TreeTools SplitsInBinaryTree
#' @export
SplitsInBinaryTree.ClusterTable <- function (tree) NTip(tree) - 3L
