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
#' Depict the splits that are matched between two trees using a specified 
#' [Generalized Robinson&ndash;Foulds](https://ms609.github.io/TreeDist/articles/Generalized-RF.html)
#' similarity measure.
#' 
#' Note that when visualizing a Robinson&ndash;Foulds distance (using 
#' `Func = RobinsonFouldsMatching`), matched splits are assigned a _similarity_
#' score of 1, which is deducted from the total number of splits to calculate 
#' the Robinson&ndash;Foulds _distance_.  Unmatched splits thus contribute one to 
#' total tree distance.
#' 
#' @param Func Function used to construct tree similarity.
#' @param tree1,tree2 Trees of class `phylo`, with identical leaf labels.
#' @param setPar Logical specifying whether graphical parameters should be 
#' set to display trees side by side.
#' @param precision Integer specifying number of significant figures to display
#' when reporting matching scores.
#' @param Plot Function to use to plot trees.
#' @param matchZeros Logical specifying whether to pair splits with zero
#' similarity (`TRUE`), or leave them unpaired (`FALSE`).
#' @param plainEdges Logical specifying whether to plot edges with a uniform
#' width and colour (`TRUE`), or whether to draw edge widths according to the
#' similarity of the associated splits (`FALSE`).
#' @param edge.width,edge.color,\dots Additional parameters to send to `Plot()`.
#' 
#' @importFrom ape nodelabels edgelabels plot.phylo
#' @importFrom colorspace qualitative_hcl sequential_hcl
#' @importFrom graphics par
#' @importFrom TreeTools as.Splits
#' 
#' @examples 
#' tree1 <- TreeTools::BalancedTree(6)
#' tree2 <- TreeTools::PectinateTree(6)
#' 
#' VisualizeMatching(RobinsonFouldsMatching, tree1, tree2)
#' VisualizeMatching(SharedPhylogeneticInfo, tree1, tree2, matchZeros = FALSE)
#' @template MRS
#' @encoding UTF-8
#' @export
VisualizeMatching <- function(Func, tree1, tree2, setPar = TRUE,
                              precision = 3L, Plot = plot.phylo,
                              matchZeros = TRUE, plainEdges = FALSE,
                              edge.width = 1, edge.color = "black",
                              ...) {
  
  splits1 <- as.Splits(tree1)
  edge1 <- tree1$edge
  child1 <- edge1[, 2]
  nTip <- attr(splits1, "nTip")
  splitEdges1 <- vapply(as.integer(rownames(splits1)),
                        function(node) which(child1 == node), integer(1))
  
  splits2 <- as.Splits(tree2, tipLabels = tree1)
  edge2 <- tree2$edge
  child2 <- edge2[, 2]
  splitEdges2 <- vapply(as.integer(rownames(splits2)),
                        function(node) which(child2 == node), integer(1))
  
  matching <- Func(tree1, tree2, reportMatching = TRUE)
  pairings <- attr(matching, "matching")
  scores <- attr(matching, "pairScores")
  pairScores <- signif(mapply(function(i, j) scores[i, j],
                              seq_along(pairings), pairings), precision)
  
  adjNo <- c(0.5, -0.2)
  adjVal <- c(0.5, 1.1)
  faint <- "#aaaaaa"
  
  if (setPar) {
    origPar <- par(mfrow = c(1, 2), mar = rep(0.5, 4))
    on.exit(par(origPar))
  }
  
  LabelUnpaired <- function(splitEdges, unpaired) {
    if (any(unpaired)) {
      #edgelabels(text="\u2012", edge=splitEdges[unpaired],
      edgelabels(text = expression("-"), edge = splitEdges[unpaired],
                 frame = "n", col = faint, adj = adjNo)
      edgelabels(text = "0", edge = splitEdges[unpaired],
                 frame = "n", col = faint, cex = 0.8, adj = adjVal)
    }
  }
  
  if (plainEdges) {
    Plot(tree1, edge.width = edge.width, edge.color = edge.color, ...)
  } else {
    Normalize <- function(scores, na.rm = FALSE) {
      if (length(scores) == 0) return(scores)
      if (na.rm) {
        scores <- scores[!is.na(scores)]
      } else {
        scores[is.na(scores)] <- 0
      }
      if (any(scores < 0)) {
        stop("Negative scores not supported")                                   # nocov
      }
      if (max(scores) == 0) return(scores)
      if (min(scores) == max(scores)) return(rep(1L, length(scores)))
      
      scores / max(scores)
    }
    
    OtherRootEdge <- function(splitNodes, edge) {
      parent <- edge[, 1]
      child <- edge[, 2]
      rootEdges <- which(parent == min(parent))
      rootChildren <- child[rootEdges]
      splitEdges <- vapply(splitNodes, match, table = child, 0)
      got <- rootChildren %in% splitNodes
      if (any(got)) {
        c(score = as.integer(which(splitNodes == rootChildren[got])),
          edge = rootEdges[!got])
      } else {
        c(score = NA, edge = NA)
      }
    }
    edgeColPalette <- sequential_hcl(n = 256L, palette = "Viridis")
     
    EdgyPlot <- function(tree, splits, edge, splitEdges,
                         normalizedScores, ...) {
      splitNodes <- as.integer(names(splits))
      ore <- OtherRootEdge(splitNodes, edge)
      if (length(normalizedScores) && !is.na(ore[1])) {
        ns <- c(normalizedScores, normalizedScores[ore["score"]])
        se <- c(splitEdges, ore[2])
      } else {
        ns <- normalizedScores
        se <- splitEdges
      }
      
      edge.width <- rep(1, nrow(edge))
      edge.width[se] <-  1 + (10 * ns)
      edge.color <- rep("black", nrow(edge))
      edge.color[se] <- edgeColPalette[1 + ceiling(255 * ns)]
      
      Plot(tree, edge.width = edge.width, edge.color = edge.color, ...)
    }
    
    EdgyPlot(tree1, splits1, edge1, splitEdges1, Normalize(pairScores), ...)
  }
  paired1 <- if (matchZeros) {
    !is.na(pairings)
  } else {
    !is.na(pairings) & pairScores > 0
  }
  palette <- qualitative_hcl(sum(paired1), c = 42, l = 88)
  pairedPairScores <- pairScores[paired1]
  pairLabels <- seq_len(sum(paired1))
  if (any(pairLabels)) {
    edgelabels(text = pairLabels, edge = splitEdges1[paired1],
               bg = palette, adj = adjNo)
    edgelabels(text = pairedPairScores, edge = splitEdges1[paired1], 
               frame = "n", adj = adjVal, cex = 0.8,
               col = ifelse(pairedPairScores, "black", faint))
  }
  LabelUnpaired(splitEdges1, !paired1)
  
  
  paired2 <- seq_along(splitEdges2) %in% pairings[paired1]
  pairNames2 <- pairings[paired1]
  
  if (plainEdges) {
    Plot(tree2, edge.width = edge.width, edge.color = edge.color, ...)
  } else {
    EdgyPlot(tree2, splits2[[pairNames2]], edge2, splitEdges2[pairNames2],
             Normalize(pairedPairScores, na.rm = TRUE), ...)
  }
  if (any(pairLabels)) {
    edgelabels(text = pairLabels, edge = splitEdges2[pairNames2],
               bg = palette, adj=adjNo)
    edgelabels(text = pairedPairScores, edge = splitEdges2[pairNames2], 
               frame = "n", adj = adjVal, cex = 0.8,
               col = ifelse(pairedPairScores, "black", faint))
  }
  LabelUnpaired(splitEdges2, !paired2)
  
  # Return:
  invisible()
}

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
  x <- mapping[, 1]
  y <- mapping[, 2]
  logStrain <- apply(mstEnds, 1, function(ends) {
    orig <- distMat[ends[1], ends[2]]
    mapped <- sum((x[ends[1]] - y[ends[2]]) ^ 2)
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
