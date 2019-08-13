#' Plot a simple tree
#'
#' Convenience plotting function used in vignettes and documentation.
#'
#' @param tr A tree of class `phylo`, with tips labelled as integers
#' @param title `main` title for the plot
#' @param bold Integer specifying which tips to print in bold
#' @param leaveRoom Logical specifying whether to leave space to print
#' tree distances beneath the plotted tree
#' @param prune,graft Integer vectors specifying which edges to highlight as
#' pruned and grafted.
#' @param \dots Additional parameters to `plot.phylo`
#'
#' @author Martin R. Smith
#' @keywords internal
#' @export
TreeDistPlot <- function (tr, title=NULL, bold=NULL, leaveRoom = TRUE,
                     prune=integer(0), graft=integer(0), ...) {
  if (is.null(tr$edge.length)) tr$edge.length <- rep(1, dim(tr$edge)[1])
  if (is.null(tr$edge.width)) tr$edge.width <- rep(1, dim(tr$edge)[1])
  edge.color <- rep('black', dim(tr$edge)[1])

  if (length(prune) > 0 || length(graft) > 0) {
    tr$edge.width[c(prune, graft)] <- 3L
    edge.color[prune] <- cbPalette8[7]
    edge.color[graft] <- cbPalette8[4]
  }

  tipNumbers <- tr$tip.label
  font <- rep(1L, length(tipNumbers))
  if (!is.null(bold)) font[tipNumbers %in% bold] <- 4L
  yLim <- if (leaveRoom) c(-0.4 - length(legendSequence), 11) else c(-0.4, 11)
  tr$tip.label <- LETTERS[as.integer(tipNumbers)]
  plot(tr, tip.col=Ternary::cbPalette15[-c(4, 7)][as.integer(tipNumbers)],
       main=title, cex.main=0.8, font=font,
       edge.width=tr$edge.width, edge.color=edge.color,
       y.lim=yLim, ...)
}
