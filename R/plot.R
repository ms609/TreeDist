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
#' @importFrom ape plot.phylo
#' @keywords internal
#' @export
TreeDistPlot <- function (tr, title=NULL, bold=NULL, leaveRoom = TRUE,
                     prune=integer(0), graft=integer(0), ...) {
  
  if (is.null(tr$edge.length)) tr$edge.length <- rep(1, dim(tr$edge)[1])
  if (is.null(tr$edge.width)) tr$edge.width <- rep(1, dim(tr$edge)[1])
  edge.color <- rep('black', dim(tr$edge)[1])

  if (length(prune) > 0 || length(graft) > 0) {
    tr$edge.width[c(prune, graft)] <- 3L
    edge.color[prune] <- "#D55E00" # Ternary::cbPalette8[7]
    edge.color[graft] <- "#009E73" # Ternary::cbPalette8[4]
  }
  tipCols <- c("#000000", "#004949", "#009292", "#FFB6DB", "#490092", "#B66DFF", 
               "#6DB6FF", "#B6DBFF", "#920000", "#924900", "#DB6D00", "#24FF24", 
               "#FFFF6D") # Ternary::cbPalette15[-c(4, 7)]

  tipNumbers <- tr$tip.label
  font <- rep(1L, length(tipNumbers))
  if (!is.null(bold)) font[tipNumbers %in% bold] <- 4L
  yLim <- if (leaveRoom) c(-0.4 - 8 # = -0.4 - length(legendSequence)
                           , 11) else c(-0.4, 11)
  tr$tip.label <- LETTERS[as.integer(tipNumbers)]
  plot.phylo(tr, tip.color=tipCols[as.integer(tipNumbers)], 
             main=title, cex.main=0.8, font=font, 
             edge.width=tr$edge.width, edge.color=edge.color, 
             y.lim=yLim, ...)
}

#' Visualise a matching
#' 
#' Depicts the bipartitions that are matched between two trees using a 
#' specified Generalized Robinson Foulds tree distance measure.
#' 
#' @param Func Function used to construct tree similarity.
#' @param tree1,tree2 Trees of class `phylo`, with tips labelled identically.
#' @param setPar Logical specifying whether graphical parameters should be 
#' set to display trees side by side.
#' @param precision Integer specifying number of significant figures to display
#' when reporting matching scores.
#' @param Plot Function to use to plot trees.
#' @param \dots Additional parameters to send to `Plot`.
#' 
#' @author Martin R. Smith
#' @importFrom ape nodelabels edgelabels plot.phylo
#' @importFrom colorspace qualitative_hcl
#' @importFrom graphics par
#' @export
VisualizeMatching <- function(Func, tree1, tree2, setPar = TRUE, 
                              precision=3, Plot = plot.phylo, ...) {
  
  splits1 <- Tree2Splits(tree1)
  edge1 <- tree1$edge
  child1 <- edge1[, 2]
  partitionEdges1 <- vapply(colnames(splits1), 
                            function (node) which(child1 == node), integer(1))
  
  splits2 <- Tree2Splits(tree2)
  edge2 <- tree2$edge
  child2 <- edge2[, 2]
  partitionEdges2 <- vapply(colnames(splits2), 
                            function (node) which(child2 == node), integer(1))
  
  matching <- Func(tree1, tree2, reportMatching = TRUE)
  pairings <- attr(matching, 'matching')
  scores <- attr(matching, 'pairScores')
  pairScores <- signif(mapply(function (i, j) scores[i, j],
                              seq_along(pairings), pairings), precision)
  
  nSplits <- length(pairings)
  palette <- qualitative_hcl(nSplits, c=42, l=88)
  adjNo <- c(0.5, -0.2)
  adjVal <- c(0.5, 1.1)
  
  if (setPar) origPar <- par(mfrow=c(2, 1), mar=rep(0.5, 4))
  
  Plot(tree1, ...)
  edgelabels(text=seq_along(pairings), edge=partitionEdges1, bg=palette, adj=adjNo)
  edgelabels(text=pairScores, edge=partitionEdges1, 
             frame='n', adj=adjVal, cex=0.8)
  
  Plot(tree2, ...)
  edgelabels(text=order(pairings), edge=partitionEdges2, bg=palette[order(pairings)], adj=adjNo)
  edgelabels(text=pairScores[order(pairings)], edge=partitionEdges2, 
             frame='n', adj=adjVal, cex=0.8)
  
  if (setPar) par(origPar)
}
