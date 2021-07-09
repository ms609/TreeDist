#' Convert phylogenetic tree to `ClusterTable`
#'
#' `as.ClusterTable()` converts a phylogenetic tree to a `ClusterTable` object,
#' which is an internal representation of its splits suitable for rapid tree
#' distance calculation (per Day, 1985).
#' 
#' Each row of a cluster table relates to a clade on a tree rooted on tip 1.
#' Tips are numbered according to the order in which they are visited in 
#' preorder: i.e., if plotted using `plot(x)`, from the top of the page
#' downwards.  A clade containing the tips 2 .. 5 would be denoted by the
#' entry `2, 5`, in either row 2 or row 5 of the cluster table.
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
#' @seealso [S3 methods][ClusterTable-methods] for `ClusterTable` objects.
#' @examples
#' tree1 <- ape::read.tree(text = "(A, (B, (C, (D, E))));");
#' tree2 <- ape::read.tree(text = "(A, (B, (D, (C, E))));");
#' ct1 <- as.ClusterTable(tree1)
#' summary(ct1)
#' as.matrix(ct1)
#' 
#' # Tip label order must match ct1 to allow comparison
#' ct2 <- as.ClusterTable(tree2, tipLabels = LETTERS[1:5])
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
as.ClusterTable.list <- function (x, tipLabels = NULL, ...) {
  lapply(x, as.ClusterTable,
         tipLabels = if (is.null(tipLabels)) TipLabels(x) else tipLabels, ...)
}

#' @rdname ClusterTable
#' @export
as.ClusterTable.multiPhylo <- as.ClusterTable.list

#' S3 methods for `ClusterTable` objects
#' 
#' S3 methods for [`ClusterTable`] objects.
#' 
#' @param x,object Object of class `ClusterTable`.
#' @param \dots Additional arguments for consistency with S3 methods.
#'
#' @examples
#' clustab <- as.ClusterTable(TreeTools::BalancedTree(6))
#' as.matrix(clustab)
#' @template MRS
#' @name ClusterTable-methods
#' @export
as.matrix.ClusterTable <- function (x, ...) {
  ClusterTable_matrix(x)
}

#' @rdname ClusterTable-methods
#' @examples
#' print(clustab)
#' @export
print.ClusterTable <- function (x, ...) {
  nTip <- attr(x, 'nTip')
  labels <- attr(x, 'tip.label')
  cat("ClusterTable on" , nTip, "leaves:", labels[1], "..", labels[nTip])
}

#' @rdname ClusterTable-methods
#' @examples
#' summary(clustab)
#' @export
summary.ClusterTable <- function (object, ...) {
  nTip <- attr(object, 'nTip')
  mat <- ClusterTable_matrix(object)
  cat("ClusterTable on" , nTip, "leaves:\n")
  cat(" ", rep(c(1:9, ' '), length.out = nTip), "\n", sep = '')
  apply(mat, 1, function (x) {
    if (x[1] > 0) {
      cat(' ', rep('.', x[1] - 1), rep('*', 1 + x[2] - x[1]),
          rep('.', nTip - x[2]), "\n", sep = '')
    }
  })
  
  cat(paste0(" ", seq_len(nTip), ": ", 
             attr(object, 'tip.label')[ClusterTable_decode(object)]), "\n")
}
