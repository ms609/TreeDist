#' Plot a simple tree
#'
#' Convenience plotting function used in vignettes and documentation.
#'
#' @param tr A tree of class `phylo`, with tips labelled as integers
#' @param title `main` title for the plot
#' @param bold Integer specifying which tips to print in bold
#' @param leaveRoom Logical specifying whether to leave space to print
#' tree distances beneath the plotted tree
#' @param edge.color Additional parameter to `plot.phylo`; will be overridden
#' by `prune` and `graft` as requested.
#' @param prune,graft Integer vectors specifying which edges to highlight as
#' pruned and grafted.
#' @param edge.width,\dots Additional parameters to `plot.phylo`
#'
#' @author Martin R. Smith
#' @importFrom ape plot.phylo
#' @keywords internal
#' @export
TreeDistPlot <- function (tr, title=NULL, bold=NULL, leaveRoom = FALSE,
                     prune=integer(0), graft=integer(0), edge.color = 'black',
                     edge.width = NULL,
                     ...) {
  
  if (is.null(tr$edge.length)) tr$edge.length <- rep(1, dim(tr$edge)[1])
  if (is.null(edge.width)) {
    edge.width <- if (is.null(tr$edge.width)) {
      rep(1, dim(tr$edge)[1])
    } else {
      tr$edge.width
    }
  }
  if (length(edge.color) == 1) edge.color <- rep(edge.color, dim(tr$edge)[1])
  nTip <- length(tr$tip.label)
  if (all(tr$tip.label %in% LETTERS)) {
    tr$tip.label <- match(tr$tip.label, LETTERS)
  } else if (all(tr$tip.label %in% letters)) {
    tr$tip.label <- match(tr$tip.label, letters)
  }

  if (length(prune) > 0 || length(graft) > 0) {
    edge.width[c(prune, graft)] <- pmax(3L, edge.width[c(prune, graft)])
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
                           , nTip) else c(-0.4, nTip)
  tr$tip.label <- LETTERS[as.integer(tipNumbers)]
  plot.phylo(tr, tip.color = tipCols[as.integer(tipNumbers)], 
             main = title, cex.main = 0.8, font = font, 
             edge.width = edge.width, edge.color = edge.color, 
             y.lim=yLim, ...)
}

#' Visualise a matching
#' 
#' Depicts the bipartitions that are matched between two trees using a 
#' specified Generalized Robinson-Foulds similarity measure.
#' 
#' Note that when visualizing a Robinson-Foulds distance (using 
#' `Func = RobinsonFoulds`), matched splits are assigned a _similarity_ score
#' of 1, which is deducted from the total number of splits to calculate the
#' Robinson-Foulds _distance_.  Unmatched splits thus contribute one to total
#' tree distance.
#' 
#' @param Func Function used to construct tree similarity.
#' @param tree1,tree2 Trees of class `phylo`, with tips labelled identically.
#' @param setPar Logical specifying whether graphical parameters should be 
#' set to display trees side by side.
#' @param precision Integer specifying number of significant figures to display
#' when reporting matching scores.
#' @param Plot Function to use to plot trees.
#' @param matchZeros Logical specifying whether to pair partitions with zero
#' similarity (`TRUE`), or leave them unpaired (`FALSE`).
#' @param plainEdges Logical specifying whether to plot edges with a uniform
#' width and colour (`TRUE`), or whether to draw edge widths according to the
#' similarity of the associated splits (`FALSE`).
#' @param edge.width,edge.color,\dots Additional parameters to send to `Plot`.
#' 
#' @author Martin R. Smith
#' @importFrom ape nodelabels edgelabels plot.phylo
#' @importFrom colorspace qualitative_hcl sequential_hcl
#' @importFrom graphics par
#' @importFrom TreeTools as.Splits
#' @export
VisualizeMatching <- function(Func, tree1, tree2, setPar = TRUE,
                              precision = 3L, Plot = plot.phylo,
                              matchZeros = TRUE, plainEdges = FALSE,
                              edge.width = NULL, edge.color = NULL,
                              ...) {
  
  splits1 <- as.logical(as.Splits(tree1))
  edge1 <- tree1$edge
  child1 <- edge1[, 2]
  partitionEdges1 <- vapply(nTip(splits2) + 1L + seq_len(splits2), 
                            function (node) which(child1 == node), integer(1))
  
  splits2 <- as.logical(as.Splits(tree2))
  edge2 <- tree2$edge
  child2 <- edge2[, 2]
  partitionEdges2 <- vapply(nTip(splits2) + 1L + seq_len(splits2), 
                            function (node) which(child2 == node), integer(1))
  
  matching <- Func(tree1, tree2, reportMatching = TRUE)
  pairings <- attr(matching, 'matching')
  scores <- attr(matching, 'pairScores')
  pairScores <- signif(mapply(function (i, j) scores[i, j],
                              seq_along(pairings), pairings), precision)
  
  adjNo <- c(0.5, -0.2)
  adjVal <- c(0.5, 1.1)
  faint <- '#aaaaaa'
  
  if (setPar) origPar <- par(mfrow=c(1, 2), mar=rep(0.5, 4))
  
  LabelUnpaired <- function (partitionEdges, unpaired) {
    if (any(unpaired)) {
      edgelabels(text='\u2212', edge=partitionEdges[unpaired],
                 frame='n', col=faint, adj=adjNo)
      edgelabels(text='0', edge=partitionEdges[unpaired],
                 frame='n', col=faint, cex=0.8, adj=adjVal)
    }
  }
  
  if (plainEdges) {
    Plot(tree1, edge.width = edge.width, edge.color = edge.color, ...)
  } else {
    Normalize <- function (scores, na.rm = FALSE) {
      if (length(scores) == 0) return(scores)
      if (na.rm) {
        scores <- scores[!is.na(scores)]
      } else {
        scores[is.na(scores)] <- 0
      }
      if (any(scores < 0)) stop ('Negative scores not supported')
      if (max(scores) == 0) return (scores)
      if (min(scores) == max(scores)) return (rep(1L, length(scores)))
      
      scores / max(scores)
    }
    
    OtherRootEdge <- function (splitNames, edge) {
      parent <- edge[, 1]
      rootEdges <- which(parent == min(parent))
      got <- edge[rootEdges, 2L] %in% splitNames
      if (any(got)) {
        
      c(score = which(splitNames == edge[rootEdges[got], 2L]),
        edge = rootEdges[!got])
      } else {
        c(score = NA, edge = NA)
      }
    }
    edgeColPalette <- sequential_hcl(n = 256L, palette = 'Viridis')
     
    EdgyPlot <- function (tree, splits, edge, partitionEdges, 
                          normalizedScores, ...) {
      ore <- OtherRootEdge(colnames(splits), edge)
      if (all(!is.na(ore)) && length(normalizedScores)) {
        ns <- c(normalizedScores, normalizedScores[ore['score']])
        pe <- c(partitionEdges, ore[2])
      } else {
        ns <- normalizedScores
        pe <- partitionEdges
      }
      
      edge.width <- rep(1, length(edge))
      edge.width[pe] <-  1 + (10 * ns)
      edge.color <- rep('black', length(edge))
      edge.color[pe] <- edgeColPalette[1 + ceiling(255 * ns)]
      
      Plot(tree, edge.width = edge.width, edge.color = edge.color, ...)
    }
    
    EdgyPlot(tree1, splits1, edge1, partitionEdges1, Normalize(pairScores), ...)
  }
  paired1 <- if (matchZeros) {
    !is.na(pairings)
  } else {
    !is.na(pairings) & pairScores > 0
  }
  palette <- qualitative_hcl(sum(paired1), c=42, l=88)
  pairedPairScores <- pairScores[paired1]
  pairLabels <- seq_len(sum(paired1))
  if (any(pairLabels)) {
    edgelabels(text=pairLabels, edge=partitionEdges1[paired1],
               bg=palette, adj=adjNo)
    edgelabels(text=pairedPairScores, edge=partitionEdges1[paired1], 
               frame='n', adj=adjVal, cex=0.8,
               col=ifelse(pairedPairScores, 'black', faint))
  }
  LabelUnpaired(partitionEdges1, !paired1)
  
  
  paired2 <- seq_along(partitionEdges2) %in% pairings[paired1]
  pairNames2 <- pairings[paired1]
  
  if (plainEdges) {
    Plot(tree2, edge.width = edge.width, edge.color = edge.color, ...)
  } else {
    EdgyPlot(tree2, splits2, edge2, partitionEdges2[pairNames2],
             Normalize(pairedPairScores, na.rm = TRUE), ...)
  }
  if (any(pairLabels)) {
    edgelabels(text=pairLabels, edge=partitionEdges2[pairNames2],
               bg=palette, adj=adjNo)
    edgelabels(text=pairedPairScores, edge=partitionEdges2[pairNames2], 
               frame='n', adj=adjVal, cex=0.8,
               col=ifelse(pairedPairScores, 'black', faint))
  }
  LabelUnpaired(partitionEdges2, !paired2)
  
  if (setPar) par(origPar)
}
