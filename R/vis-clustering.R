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


#' Visualize clustering structure via Multi-dimensional scaling
#' 
#' @param maxIter Integer specifying maximum number of moves to propose.
#' @param stability After proposing `stability` moves without improvement,
#' stop.
#' 
#' @examples
#' library("TreeTools", quietly = TRUE)
#' treeNumbers <- c(0:20, 40:60, 120:140, 200:220)
#' trees <- as.phylo(treeNumbers, 8)
#' distances <- ClusteringInfoDist(trees)
#' clustering <- kmeans(distances, 4)
#' clust <- clustering$cluster
#' cmds <- cmdscale(distances)
#' 
#' par(mfrow = c(2, 3))
#' plot(cmds, col = clust, cex = 1.5, asp = 1, axes = FALSE, ann = FALSE)
#' plot(cluster::silhouette(clust, distances), col = 1:4, main = "Original")
#' plot(cluster::silhouette(clust, dist(cmds)), col = 1:4, main = "Classic MDS")
#' 
#' silmds <- ClusterMDS(clust, distances, cmds)
#' plot(silmds, col = clust, cex = 1.5, asp = 1, axes = FALSE, ann = FALSE)
#' MSTSegments(silmds, MSTEdges(distances),
#'             col = StrainCol(distances, silmds))
#' plot(cluster::silhouette(clust, dist(silmds)), col = 1:4, main = "SilMDS")
#' 
#' 
#' @template MRS
#' @importFrom cli cli_progress_bar cli_progress_done cli_progress_update
#' @importFrom TreeTools MSTEdges
#' @export
ClusterMDS <- function(clust, distances, proposal = cmdscale(distances),
                       maxIter = 1e6, stability = 42,
                       weight = c("sil" = 1,
                                  "mst" = 1,
                                  "sor" = 1,
                                  "sov" = 1,
                                  "mcd" = 1,
                                  "mnn" = 1,
                                  "mse" = 1,
                                  "icd" = 0),
                       trace = FALSE) {
  distances <- as.dist(distances)
  treeID <- combn(attr(distances, "Size"), 2)
  Measuring <- function(x) !is.na(weight[x]) && weight[x] != 0
  Score <- function(orig, new) sqrt(sum((orig - new) ^ 2))
  scores <- 0 * weight
  
  if (Measuring("icd")) {
    clustID <- matrix(clust[treeID], 2)
    nClust <- max(clust)
    clustPairs <- combn(nClust + 1, 2)
    clustPairs[2, ] <- clustPairs[2, ] - 1L
    clustPairKeys <- apply(clustPairs, 2, paste, collapse = "-")
    clustKey <- factor(apply(clustID, 2, 
                             function(x) paste(sort(x), collapse = "-")),
                       clustPairKeys)
    clustPair <- vapply(seq_along(clustPairKeys), function(pair) {
      clustKey == clustPairKeys[pair]
    }, logical(length(clustKey)))
    interClustDists <- apply(clustPair, 2, function(x) mean(distances[x]))
  }
  
  if (Measuring("mst")) {
    mst <- MSTEdges(distances)
    onMST <- duplicated(cbind(t(mst), treeID), MARGIN = 2)[
      -seq_len(dim(mst)[1])]
    mstLength <- distances[onMST]
    mstLength <- mstLength / sum(mstLength)
  }
  
  if (Measuring("mnn")) {
    mnnZero <- MeanNN(distances)
  }
  
  sil <- cluster::silhouette(clust, distances)
  silWidth <- sil[, "sil_width"]
  bestScore <- Inf
  mover <- 0L
  accepted <- logical(maxIter)
  la <- integer(maxIter)
  d2Log <- double(maxIter)
  icdLog <- double(maxIter)
  mstLog <- double(maxIter)
  lastAccepted <- 0L
  
  
  range <- apply(proposal, 2, range)
  scale <- range[2, ] - range[1, ]
  
  cli_progress_bar("Rescaling", total = maxIter, .auto_close = FALSE)
  for (i in seq_len(maxIter)) {
    dprime <- dist(proposal)
    
    if (Measuring("icd")) {
      icd <- apply(clustPair, 2, function(x) sum(dprime[x]) / sum(x))
      rescale <- sum(interClustDists) / sum(icd)
      proposal <- proposal * rescale
      icd <- icd * rescale
      scores["icd"] <- sqrt(sum((interClustDists - icd) ^ 2))
    }
    
    if (Measuring("mst")) {
      mstPrime <- dprime[onMST]
      mstPrime <- mstPrime / sum(mstPrime)
      scores["mst"] <- Score(mstLength, mstPrime)
    }
    
    if (Measuring("mnn")) {
      scores["mnn"] <- Score(mnnZero, MeanNN(proposal, cluster))
    }


    mappedSil <- cluster::silhouette(clust, dprime)
    dSil <- mappedSil[, "sil_width"] - silWidth
    dSil2 <- dSil * dSil
    d2Cum <- cumsum(dSil2)
    d2Sum <- d2Cum[length(d2Cum)]
    
    score <- (d2Sum * weight[["sil"]]) + icdScore + mstScore
    
    if (score < bestScore) {
      cli_progress_update(set = i)#,
                          #extra = paste0("Stability: ", lastAccepted))
      mapped <- proposal
      bestScore <- score
      accepted[i] <- TRUE
      lastAccepted <- 0L
    }
    
    d2Log[i] <- d2Sum
    icdLog[i] <- icdScore
    mstLog[i] <- mstScore
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
  
  if (trace) {
    plot(cumsum(accepted[seq_len(i)]), type = 'n')
    # points(la[seq_len(i)] * sum(accepted) / stability, type = 'l', col = 3)
    points(d2Log[seq_len(i)] * sum(accepted) / max(d2Log),
           col = 4, type = "l")
    points(mstLog[seq_len(i)] * sum(accepted) / max(mstLog),
           col = 5, type = "l")
    points(icdLog[seq_len(i)] * sum(accepted) / max(icdLog),
           col = 6, type = "l")
    points(cumsum(accepted[seq_len(i)]), type = 'l', col = 2)
    legend("top", c("Accepted", "lastAccept", "sil-d2", "mst", "icd"),
           col = 2:6, pch = 15, bty = "n", pt.cex = 2.5)
  }
  mapped
}
