#' Add minimum spanning tree to plot, colouring by stress
#' 
#' To identify strain in a multidimensional scaling of distances, it can be
#' useful to plot a minimum spanning tree
#' \insertCite{Gower1966,SmithSpace}{TreeDist}.  Colouring each edge of the
#' tree according to its strain can identify areas where the mapping is
#' stretched or compressed.
#' 
#' @param mapping Two-column matrix giving _x_ and _y_ coordinates of plotted
#' points.
#' @param mstEnds Two-column matrix identifying rows of `mapping` at end of
#' each edge of the MST, as output by [`TreeTools::MSTEdges()`].
#' @param distances Matrix or `dist` object giving original distances between
#' each pair of points.
#' @param palette Vector of colours with which to colour edges.
#' @param \dots Additional arguments to [`segments()`].
#'
#' @examples
#' set.seed(0)
#' library("TreeTools", quietly = TRUE)
#' distances <- ClusteringInfoDist(as.phylo(5:16, 8))
#' mapping <- cmdscale(distances, k = 2)
#' mstEnds <- MSTEdges(distances)
#' 
#' # Set up blank plot
#' plot(mapping, asp = 1, frame.plot = FALSE, ann = FALSE, axes = FALSE,
#'      type = "n")
#' # Add MST
#' MSTSegments(mapping, mstEnds, 
#'             col = StrainCol(distances, mapping, mstEnds))
#' # Add points at end so they overprint the MST
#' points(mapping)
#' PlotTools::SpectrumLegend(
#'  "bottomleft",
#'  legend = c("Extended", "Median", "Contracted"),
#'  bty = "n",     # No box
#'  y.intersp = 2, # Expand in Y direction
#'  palette = hcl.colors(256L, "RdYlBu", rev = TRUE)
#' )
#' @template MRS
#' @references \insertAllCited{}
#' @family tree space functions
#' @importFrom graphics segments
#' @export
MSTSegments <- function(mapping, mstEnds, ...) {
  segments(mapping[mstEnds[, 1], 1], mapping[mstEnds[, 1], 2],
           mapping[mstEnds[, 2], 1], mapping[mstEnds[, 2], 2], ...)
}

#' @rdname MSTSegments
#' @return `StrainCol()` returns a vector in which each entry is selected from
#' `palette`, with an attribute `logStrain` denoting the logarithm of the
#' mapped over original distance, shifted such that the median value is zero.
#' Palette colours are assigned centred on the median value, with entries
#' early in `palette` assigned to edges in which the ratio of mapped
#' distance to original distance is small.
#' @importFrom grDevices hcl.colors
#' @importFrom TreeTools MSTEdges
#' @export
StrainCol <- function(distances, mapping, mstEnds = MSTEdges(distances),
                      palette = rev(hcl.colors(256L, "RdYlBu"))) {
  distMat <- as.matrix(distances)
  logStrain <- apply(mstEnds, 1, function(ends) {
    orig <- distMat[ends[1], ends[2]]
    mapped <- sum((mapping[ends[1], ] - mapping[ends[2], ]) ^ 2)
    (
      log(mapped) / 2 # sqrt
    ) - log(orig) # High when mapping extends original distances
  })
  strain <- logStrain - median(logStrain[is.finite(logStrain)])
  # Infinite values arise when orig == 0
  maxVal <- max(abs(strain[is.finite(strain)])) + sqrt(.Machine$double.eps)
  nCols <- length(palette)
  bins <- cut(strain, seq(-maxVal, maxVal, length.out = nCols))
  
  # Return:
  structure(palette[bins],
            logStrain = strain)
}
