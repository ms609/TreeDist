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
#' 
#' # Tip label order must match ct1 to allow comparison
#' ct2 <- as.ClusterTable(tree2, tipLabels = LETTERS[1:5])
#' 
#' #TODO: RF (ct1, ct2)
#' @template MRS
#'
#' @export
as.ClusterTable <- function (x, tipLabels = NULL, ...) UseMethod('as.ClusterTable')

#' @rdname as.ClusterTable
#' @importFrom TreeTools NTip RenumberTips
#' @export
as.ClusterTable.phylo <- function (x, tipLabels = NULL, ...) {
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

#' Print `ClusterTable` object
#'
#' S3 method for objects of class `ClusterTable`.
#'
#' @param x Object of class `ClusterTable`.
#' @param \dots Additional arguments for consistency with S3 method (unused).
#'
#' @export
print.ClusterTable <- function (x, ...) {
  cat("ClusterTable on" , paste0(attr(x, 'nTip')), "leaves:", paste0(attr(x, 'tip.label')))
}