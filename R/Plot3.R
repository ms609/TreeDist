#' Pseudo-3D plotting
#' 
#' `Plot3()` displays three-dimensional data in two dimensions, reflecting the
#' third dimension with point scaling, overlap and fogging.
#' Points with a lower `z` value are smaller than, fainter than, and overlapped
#' by points with a higher value.
#' 
#' @param x,y,z Coordinates of points to plot.
#' @param fog Numeric from zero (no fading) to one (furthest points are
#' invisible) specifying amount to fade distant points.
#' @param shrink Numeric specifying degree to which size of plotted point
#' should reflect `z` position. `0` denotes no scaling; if `1`, furthest
#' point will have zero size.
#' @param plot.bg Colour with which to fill plot area, used as fog colour.
#' @param bg,cex,col,pch,add,axes,frame.plot,\dots Parameters passed to 
#' [plot.default()].
#' 
#' @examples 
#' Plot3(1:10, 1:10, 1:10, cex = 7, pch = 16, axes = FALSE, asp = 1)
#' # Extreme values of fog and shrink will cause smallest z values to
#' # become invisible.
#' Plot3(1:10, 1:10, 1:10, cex = 7, pch = 16, axes = FALSE, asp = 1,
#'       fog = 1, shrink = 1)
#' @template MRS
#' @importFrom graphics plot points rect
#' @importFrom grDevices colorRamp rgb
#' @export
Plot3 <- function(x, y = NULL, z = NULL,
                   pch = par("pch"), col = par("col"),
                   bg = NA, cex = 1,
                   axes = TRUE,
                   frame.plot = axes,
                   plot.bg = NA, 
                   fog = 1/2, shrink = 1/2,
                   add = FALSE,
                   ...) {
  if (is.null(y)) {
    z <- x[, 3]
    y <- x[, 2]
    x <- x[, 1]
  }
  
  zResolution <- 128L
  n <- length(z)
  zOrder <- order(z)
  bg <- rep_len(bg, n)
  col <- rep_len(col, n)
  cex <- rep_len(cex, n)
  pch <- rep_len(pch, n)
  
  #zStep <- as.integer(cut(z, zResolution))
  zScale <- max(z) - z
  zScale <- zScale / max(zScale)
  
  fogOffset <- zResolution * fog
  bgCol <- if (is.na(plot.bg)) "white" else plot.bg
  .FadeCol <- function(x, fadeAmount) {
    if (is.na(x)) {
      NA_character_
    } else {
      rgb(colorRamp(c(x, bgCol), alpha = TRUE)(fog * fadeAmount),
          maxColorValue = 255)
    }
  }
  fadedCol <- vapply(seq_along(z), function(i) {
    .FadeCol(col[i], zScale[i])
  }, character(1))
  fadedBg <- vapply(seq_along(z), function(i) {
    .FadeCol(bg[i], zScale[i])
  }, character(1))
  if (!add) {
    plot(x, y, type = "n", axes = axes, frame.plot = frame.plot, ...)
    if (!is.na(plot.bg)) {
      rect(par("usr")[[1]], par("usr")[[3]], par("usr")[[2]], par("usr")[[4]],
           col = plot.bg, border = frame.plot)
    }
  }
  points(x[zOrder], y[zOrder],
       cex = cex[zOrder] * (1 - (shrink * zScale[zOrder])),
       pch = pch[zOrder],
       col = fadedCol[zOrder],
       bg = fadedBg[zOrder],
       ...)
}
