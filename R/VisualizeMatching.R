#' Visualize a matching
#' 
#' Depict the splits that are matched between two trees using a specified 
#' [Generalized Robinson&ndash;Foulds](
#' https://ms609.github.io/TreeDist/articles/Generalized-RF.html)
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
#' @param edge.cex Character expansion for edge labels.
#' If `FALSE`, suppress edge labels.
#' @param value.cex Character expansion for values on edge labels.
#' If `FALSE`, values are not displayed.
#' @param edge.frame Character specifying the kind of frame to be printed around
#' the text of the edge labels.  Choose an abbreviation of `"rect"`, `"circle"`,
#' or `"none"`.
#' @param edge.width,edge.color,\dots Additional parameters to send to `Plot()`.
#' 
#' @returns `VisualizeMatching()` invisibly returns the matching of splits
#' between `tree1` and `tree2` (i.e. 
#' `Func(tree1, tree2, reportMatching = TRUE)`)
#' 
#' @examples 
#' tree1 <- TreeTools::BalancedTree(6)
#' tree2 <- TreeTools::PectinateTree(6)
#' 
#' VisualizeMatching(RobinsonFouldsMatching, tree1, tree2)
#' matching <- VisualizeMatching(SharedPhylogeneticInfo, tree1, tree2,
#'                               matchZeros = FALSE)
#' attributes(matching)
#' @template MRS
#' @encoding UTF-8
#' @importFrom ape nodelabels edgelabels plot.phylo
#' @importFrom colorspace qualitative_hcl sequential_hcl
#' @importFrom graphics par
#' @importFrom TreeTools as.Splits
#' @export
VisualizeMatching <- function (Func, tree1, tree2, setPar = TRUE,
                               precision = 3L, Plot = plot.phylo,
                               matchZeros = TRUE, plainEdges = FALSE,
                               edge.cex = par("cex"),
                               value.cex = edge.cex * 0.8,
                               edge.frame = "rect",
                               edge.width = 1, edge.color = "black",
                               ...) 
{
  splits1 <- as.Splits(tree1)
  edge1 <- tree1[["edge"]]
  child1 <- edge1[, 2]
  nTip <- attr(splits1, "nTip")
  splitEdges1 <- vapply(as.integer(rownames(splits1)),
                        function(node) which(child1 == node), integer(1))
  
  splits2 <- as.Splits(tree2, tipLabels = tree1)
  edge2 <- tree2[["edge"]]
  child2 <- edge2[, 2]
  splitEdges2 <- vapply(as.integer(rownames(splits2)),
                        function(node) which(child2 == node), integer(1))
  
  matching <- Func(tree1, tree2, reportMatching = TRUE)
  pairings <- attr(matching, "matching")
  scores <- attr(matching, "pairScores")
  pairScores <- signif(mapply(function(i, j) scores[i, j],
                              seq_along(pairings), pairings), precision)
  
  faint <- "#aaaaaa"
  
  if (setPar) {
    origPar <- par(mfrow = c(1, 2), mar = rep(0.5, 4))
    on.exit(par(origPar))
  }
  
  .LabelEdge <- function(label, edges, frame = "n", ...) {
    if (edge.cex) {
      edgelabels(text = label, edge = edges, frame = frame,
                 cex = edge.cex, adj = c(0.5, -0.2), ...)
    }
  }
  .LabelValue <- function(label, edges, frame = "n", ...) {
    if (value.cex) {
      edgelabels(text = label, edge = edges, frame = frame,
                 cex = value.cex, adj = c(0.5, 1.1), ...)
    }
  }
  
  .LabelUnpaired <- function(splitEdges, unpaired) {
    if (any(unpaired)) {
      .LabelEdge(label = expression("-"), edges = splitEdges[unpaired],
                 frame = "n", col = faint)
      .LabelValue(label = "0", edges = splitEdges[unpaired],
                  frame = "n", col = faint)
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
      if (max(scores) == 0) {
        return(scores)
      }
      if (min(scores) == max(scores)) {
        return(rep(1L, length(scores)))
      }
      scores / max(scores)
    }
    
    OtherRootEdge <- function(splitNodes, edge) {
      parent <- edge[, 1]
      child <- edge[, 2]
      rootEdges <- which(parent == min(parent))
      rootChildren <- child[rootEdges]
      treeIsRooted <- length(rootChildren) < 3
      if (treeIsRooted) {
        splitEdges <- vapply(splitNodes, match, table = child, 0)
        got <- rootChildren %in% splitNodes
        if (any(got)) {
          if (sum(got) != 1) {
            warning("Unexpected polytomy")
          }
          c(score = as.integer(which(splitNodes %in% rootChildren[got])),
            edge = rootEdges[!got])
        } else {
          # `edge` is not a root edge
          c(score = NA_integer_, edge = NA_integer_)
        }
      } else {
        # Tree is unrooted => there is no root edge at all
        c(score = NA_integer_, edge = NA_integer_)
      }
    }
    edgeColPalette <- sequential_hcl(n = 256L, palette = "Viridis")
    
    EdgyPlot <- function(tree, splits, edge, splitEdges,
                         normalizedScores, ...) {
      splitNodes <- as.integer(names(splits))
      ore <- OtherRootEdge(splitNodes, edge)
      if (length(normalizedScores) && !is.na(ore[[1]])) {
        ns <- c(normalizedScores, normalizedScores[ore["score"]])
        se <- c(splitEdges, ore[[2]])
      } else {
        ns <- normalizedScores
        se <- splitEdges
      }
      
      edge.width <- rep(1, nrow(edge))
      edge.width[se] <- 1 + (10 * ns)
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
    .LabelEdge(pairLabels, splitEdges1[paired1], frame = edge.frame,
               bg = palette)
    .LabelValue(pairedPairScores, splitEdges1[paired1],
                col = ifelse(pairedPairScores, "black", faint))
  }
  .LabelUnpaired(splitEdges1, !paired1)
  
  
  paired2 <- seq_along(splitEdges2) %in% pairings[paired1]
  pairNames2 <- pairings[paired1]
  
  if (plainEdges) {
    Plot(tree2, edge.width = edge.width, edge.color = edge.color, ...)
  } else {
    EdgyPlot(tree2, splits2[[pairNames2]], edge2, splitEdges2[pairNames2],
             Normalize(pairedPairScores, na.rm = TRUE), ...)
  }
  if (any(pairLabels)) {
    .LabelEdge(pairLabels, splitEdges2[pairNames2], frame = edge.frame,
               bg = palette)
    .LabelValue(pairedPairScores, splitEdges2[pairNames2],
               col = ifelse(pairedPairScores, "black", faint))
  }
  .LabelUnpaired(splitEdges2, !paired2)
  
  # Return:
  invisible(matching)
}
