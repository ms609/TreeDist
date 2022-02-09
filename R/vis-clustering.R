#' Visualize clustering distances
#' 
#' 
#' 
#' @param distances Matrix or `dist` object giving distance between each
#' pair of trees.
#' @param clust Integer vector listing the sequentially numbered cluster to
#' which each tree is assigned.
#' @param col Colours with which to colour each cluster.
#' Either specify a vector of colours with an entry for each cluster,
#' or a function that will generate a palette when passed a variable denoting
#' the number of colours required.
#' 
#' @examples
#' library("TreeTools")
#' treeNumbers <- c(0:20, 40:60, 120:140, 200:220)
#' group <- rep(1:4, each = 21)
#' trees <- as.phylo(treeNumbers, 8)
#' distances <- ClusteringInfoDist(trees)
#' set.seed(1) # Bad clustering
#' clustering <- kmeans(distances, 4)
#' plot(cmdscale(distances), col = clustering$cluster, cex = 1.5,
#'      asp = 1, axes = FALSE, ann = FALSE)
#' ClusVis(distances, clustering$cluster)
#' plot(cluster::silhouette(clustering$cluster, distances),
#'      col = 1:4)
#'
#' set.seed(2) # Good clustering
#' clustering <- kmeans(distances, 4)
#' Cols <- function (n) hcl.colors(n, pal = "Dynamic")
#' plot(cmdscale(distances), col = Cols(4)[clustering$cluster], cex = 1.5,
#'      asp = 1, axes = FALSE, ann = FALSE)
#' 
#' ClusVis(distances, clustering$cluster, col = Cols)
#' @template MRS
#' @export
#' 
ClusVis <- function(distances, clust, col = seq_len) {
  distances <- as.dist(distances)
  treeID <- combn(attr(distances, "Size"), 2)
  
  clustID <- matrix(clust[treeID], 2)
  nClust <- max(clust)
  clustPairs <- combn(nClust + 1, 2)
  clustPairs[2, ] <- clustPairs[2, ] - 1L
  clustPairKeys <- apply(clustPairs, 2, paste, collapse = " ")
  clustKey <- factor(apply(clustID, 2, 
                           function(x) paste(sort(x), collapse = " ")),
                     clustPairKeys)
  master <- hist(distances, plot = FALSE)
  breaks <- master$breaks
  counts <- vapply(seq_along(clustPairKeys), function(pair) {
    hist(distances[clustKey == clustPairKeys[pair]],
         breaks = breaks, plot = FALSE)$counts
  }, master$counts)
  clToCl <- vapply(seq_along(clustPairKeys), function(pair) {
    mean(distances[clustKey == clustPairKeys[pair]])
  }, double(1))
  
  mar <- par("mar")
  oma <- par("oma")
  oma[1] <- oma[1] + mar[1]
  mar[c(1, 3)] <- 1.1
  oPar <- par(mfrow = c(nClust, 1), mar = mar, oma = oma)
  if (is.function(col)) {
    col <- col(nClust)
  }
  on.exit(par(oPar))
  
  for (i in seq_len(nClust)) {
    thisClust <- apply(clustPairs == i, 2, any)
    clOrder <- order(clToCl[thisClust])
    countI <- counts[, thisClust][, clOrder]
    barplot(t(countI), col = col[clOrder], space = 0, ylab = "Frequency")
    legend("topleft", paste("Cluster", i),
           bty = "n", pch = 15, pt.cex = 2, col = col[i])
  }
  axis(1, at = seq_along(breaks) - 1L, breaks)
  
}


#' @examples
#' treeNumbers <- c(0:20, 40:60, 120:140, 200:220)
#' trees <- as.phylo(treeNumbers, 8)
#' distances <- ClusteringInfoDist(trees)
#' clustering <- kmeans(distances, 4)
#' clust <- clustering$cluster
#' cmds <- cmdscale(distances)
#' 
#' par(mfrow = c(2, 3))
#' plot(cmds, col = clust, cex = 1.5, asp = 1, axes = FALSE, ann = FALSE)
#' plot(cluster::silhouette(clust, distances), col = 1:4, main = "True")
#' plot(cluster::silhouette(clust, dist(cmds)), col = 1:4, main = "Classic MDS")
#' 
#' silmds <- SilMatch(clust, distances, cmds)
#' plot(silmds, col = clust, cex = 1.5, asp = 1, axes = FALSE, ann = FALSE)
#' plot(cluster::silhouette(clust, dist(silmds)), col = 1:4, main = "SilMDS")
#' 
#' 
#' @template MRS
#' @importFrom cli cli_progress_bar cli_progress_done cli_progress_update
#' @export
SilMatch <- function(clust, distances, proposal = cmdscale(distances),
                      maxIter = 1e6, stability = 1e2) {
  
  sil <- cluster::silhouette(clust, distances)
  silWidth <- sil[, "sil_width"]
  bestD2 <- Inf
  mover <- 0L
  accepted <- logical(maxIter)
  la <- integer(maxIter)
  score <- double(maxIter)
  lastAccepted <- 0L
  
  
  cli_progress_bar("Rescaling", total = maxIter, .auto_close = FALSE)
  for (i in seq_len(maxIter)) {
    dprime <- dist(proposal)
    range <- apply(proposal, 2, range)
    scale <- range[2, ] - range[1, ]
    
    mappedSil <- cluster::silhouette(clust, dprime)
    dSil <- mappedSil[, "sil_width"] - silWidth
    d2 <- dSil * dSil
    d2Cum <- cumsum(d2)
    d2Sum <- d2Cum[length(d2Cum)]
    if (d2Sum < bestD2) {
      cli_progress_update(set = i)
      mapped <- proposal
      bestD2 <- d2Sum
      accepted[i] <- TRUE
      lastAccepted <- 0L
    }
    
    score[i] <- d2Sum
    lastAccepted <- lastAccepted + 1L
    la[i] <- lastAccepted
    
    if (lastAccepted > stability) {
      break
    }
    mover <- which.max(d2Cum > runif(1, max = d2Sum))
    moveSize <- 0.2 / lastAccepted
    proposal <- mapped
    proposal[mover, ] <- proposal[mover, ] + rnorm(2, 0, moveSize * scale)
  }
  cli_progress_done()
  
  plot(cumsum(accepted[seq_len(i)]), type = 'l', col = 2)
  points(la[seq_len(i)] * sum(accepted) / stability, type = 'l', col = 3)
  points(score[seq_len(i)] * sum(accepted) / max(score), col = 4, type = "l")
  mapped
}
