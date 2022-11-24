#' Multi-dimensional scaling with distance ranges
#' 
#' @param distMin,distMax Symmetrical matrices specifying the minimum and
#' maximum possible distances between the objects to map.
#' @param maxIter Integer specifying maximum number of moves to propose.
#' @param stability After proposing `stability` moves without improvement,
#' stop.
#' 
#' @examples
#' library("TreeTools", quietly = TRUE)
#' treeNumbers <- c(0:20, 40:60, 120:140, 200:220)
#' trees <- as.phylo(treeNumbers, 8)
#' noBCD <- which(treeNumbers %in% c(10:15, 50:55, 130:135, 210:215))
#' noF <- 21 * 1:4
#' trees[noBCD] <- DropTip(trees[noBCD], c(2:4))
#' trees[noF] <- DropTip(trees[noF], 6)
#' 
#' distances <- ClusteringInfoDist(trees)
#' maxEnt <- ClusteringEntropy(trees)
#' maxObservable <- outer(maxEnt, maxEnt, "+")
#' nTip <- length(unique(unlist(TipLabels(trees))))
#' # Upper bound on maximum possible distance
#' maxPossible <- ClusteringEntropy(PectinateTree(nTip)) * 2
#' 
#' distUB <- as.matrix(distances) + maxPossible - maxObservable
#' 
#' cmds <- cmdscale(distances)
#' silmds <- ClusterMDS(clust, distances, silmds, trace = TRUE, stab = 72,
#' weight = c("sil" = 1,
#'            "ave" = 1,
#'            "max" = 1,
#'            "ks" = 1,
#'            "sor" = 1,
#'            "sov" = 1,
#'            "dfm" = 1,
#'            "mnn" = 1,
#'            "ste" = 1)
#' )
#' plot(cluster::silhouette(clust, distances), col = 1:nCl, main = "Original")
#' ClusVis(distances, clust, setPar = FALSE)
#' 
#' plot(silmds, col = clust, cex = 1.5, asp = 1, axes = FALSE, ann = FALSE)
#' plot(cluster::silhouette(clust, dist(silmds)), col = 1:nCl, main = "SilMDS")
#' ClusVis(dist(silmds), clust, setPar = FALSE)
#' 
#' plot(cmds, col = clust, cex = 1.5, asp = 1, axes = FALSE, ann = FALSE)
#' plot(cluster::silhouette(clust, dist(cmds)), col = 1:nCl, main = "Classic MDS")
#' ClusVis(dist(cmds), clust, setPar = FALSE)
#' 
#' @template MRS
#' @importFrom cli cli_progress_bar cli_progress_done cli_progress_update
#' @importFrom TreeTools MSTEdges
#' @export
MapRanges <- function(distMin, distMax, proposal = cmdscale(distMin),
                      maxIter = 1e6, stability = 42,
                      trace = FALSE) {
  lowerBound <- as.dist(distMin)
  upperBound <- as.dist(distMax)
  treeID <- combn(attr(distances, "Size"), 2)
  
  bestScore <- Inf
  mover <- 0L
  accepted <- logical(maxIter)
  la <- integer(maxIter)
  mostStable <- 0
  log <- matrix(0, maxIter, length(weight),
                dimnames = list(NULL, names(scores)))
  
  for (i in seq_len(maxIter)) {
    dprime <- dist(proposal)
    
    under <- dprime < lowerBound
    over <- dprime > upperBound
    
    score <- sqrt(sum(dprime[under | over] ^ 2))
    
    if (score < bestScore) {
      cli_progress_update(set = mostStable)
      mapped <- proposal
      bestScore <- score
      accepted[i] <- TRUE
      log[i, ] <- scores
      lastAccepted["any"] <- 0L
      lastAccepted[strategy] <- 0L
    } else {
      log[i, ] <- NA
    }
    
    lastAccepted["any"] <- lastAccepted["any"] + 1L
    lastAccepted[strategy] <- lastAccepted[strategy] + 1L
    la[i] <- lastAccepted["any"]
    mostStable <- max(mostStable, la[i])
    
    if (mostStable > stability) {
      break
    }
    
    moveSize <- 0.2 / lastAccepted[strategy]
    proposal <- mapped
    proposal[mover, ] <- t(t(proposal[mover, ]) + rnorm(2, 0, moveSize * scale))
  }
  cli_progress_done()
  
  if (trace) {
    TracePoints <- function(x) {
      if (x %in% names(scores)) {
        style <- match(x, names(scores)) + 1
        values <- log[seq_len(i), x]
        points(
          which(!is.na(values)),
          values[!is.na(values)] * sum(accepted) / max(log[, x], na.rm = TRUE),
          col = style,
          lty = (style %/% 8) + 1,
          type = "l"
        )
      }
    }
    plot(cumsum(accepted[seq_len(i)]), type = "n", frame.plot = FALSE)
    points(cumsum(accepted[seq_len(i)]), type = "l", col = 1)
    lapply(names(scores), TracePoints)
    #points(la[seq_len(i)] * sum(accepted) / stability, type = "l", col = 2)
    legend("left", c("Accepted", names(scores)),
           col = 1:15, pch = 15, bty = "n", pt.cex = 2.5, text.font = 2, cex = 1)
  }
  mapped
}
