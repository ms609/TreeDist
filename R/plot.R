#' Plot a simple tree
#'
#' Convenience plotting function used in vignettes and documentation.
#'
#' @param tr A tree of class `phylo`, with leaves labelled as integers.
#' @param title `main` title for the plot.
#' @param bold Integer specifying which leaves to print in bold.
#' @param leaveRoom Logical specifying whether to leave space to print
#' tree distances beneath the plotted tree.
#' @param edge.color Additional parameter to `plot.phylo`; will be overridden
#' by `prune` and `graft` as requested.
#' @param prune,graft Integer vectors specifying which edges to highlight as
#' pruned and grafted.
#' @param edge.width,\dots Additional parameters to `plot.phylo`.
#'
#' @template MRS
#' @importFrom ape plot.phylo
#' @keywords internal
#' @export
TreeDistPlot <- function(tr, title = NULL, bold = NULL, leaveRoom = FALSE,
                         prune = integer(0), graft = integer(0),
                         edge.color = "black", edge.width = NULL, ...) {
  
  nEdge <- dim(tr[["edge"]])[1]
  if (is.null(tr[["edge.length"]])) {
    tr[["edge.length"]] <- rep(1, nEdge)
  }
  if (is.null(edge.width)) {
    edge.width <- if (is.null(tr[["edge.width"]])) {
      rep(1, nEdge)
    } else {
      tr[["edge.width"]]
    }
  }
  if (length(edge.color) == 1) {
    edge.color <- rep(edge.color, nEdge)
  }
  nTip <- length(tr[["tip.label"]])
  if (all(tr[["tip.label"]] %in% LETTERS)) {
    tr[["tip.label"]] <- match(tr[["tip.label"]], LETTERS)
  } else if (all(tr[["tip.label"]] %in% letters)) {
    tr[["tip.label"]] <- match(tr[["tip.label"]], letters)
  }

  if (length(prune) > 0 || length(graft) > 0) {
    edge.width[c(prune, graft)] <- pmax(3L, edge.width[c(prune, graft)])
    edge.color[prune] <- "#D55E00" # Ternary::cbPalette8[7]
    edge.color[graft] <- "#009E73" # Ternary::cbPalette8[4]
  }
  tipCols <- c("#000000", "#004949", "#009292", "#FFB6DB", "#490092", "#B66DFF",
               "#6DB6FF", "#B6DBFF", "#920000", "#924900", "#DB6D00", "#24FF24",
               "#FFFF6D") # Ternary::cbPalette15[-c(4, 7)]

  tipNumbers <- tr[["tip.label"]]
  font <- rep(1L, length(tipNumbers))
  if (!is.null(bold)) {
    font[tipNumbers %in% bold] <- 4L
  }
  yLim <- if (leaveRoom) {
    c(-0.4 - 8, # = -0.4 - length(legendSequence)
      nTip) 
  } else {
    c(-0.4, nTip)
  }
  tipInts <- tryCatch(as.integer(tipNumbers),
                      warning = function(e) {
                        warning("Leaves of `tr` must be labelled with integers")
                      })
  
  tr[["tip.label"]] <- LETTERS[tipInts]
  plot.phylo(tr, tip.color = tipCols[tipInts],
             main = title, cex.main = 0.8, font = font, 
             edge.width = edge.width, edge.color = edge.color, 
             y.lim=yLim, ...)
}
