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
#' @param setPar Logical specifying whether to set graphical parameters to
#' display all plots at once.
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
ClusVis <- function(distances, clust, col = seq_len, setPar = TRUE) {
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
  
  if (setPar) {
    mar <- par("mar")
    oma <- par("oma")
    oma[1] <- oma[1] + mar[1]
    mar[c(1, 3)] <- 1.1
    oPar <- par(mfrow = c(nClust, 1), mar = mar, oma = oma)
    on.exit(par(oPar))
  }
  
  if (is.function(col)) {
    col <- col(nClust)
  }
  
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
#' @param weight A named vector listing weights to apply to clustering
#' properties when constructing mapping. Options:
#' `sil`: silhouette coefficients
#' `ave`: Mean and median distances between points within each cluster,
#' and between each pair of clusters
#' `max`: Minimum and maximum distances between points within each cluster,
#' and between each pair of clusters
#' 
#' `ks`: Kolmogorov–Smirnov test statistic between original and mapped
#' distances, conducted for all cluster pairs
#' `sor`: Sum of Ranges of each cluster
#' `sov`: Sum of Variance within each cluster
#' `dfm`: Mean distance from median
#' `mnn`: Mean distance to nearest neighbour
#' `ste`: Mean length of minimum spanning tree edge
#' 
#' @examples
#' library("TreeTools", quietly = TRUE)
#' treeNumbers <- c(0:20, 40:60, 120:140, 200:220)
#' trees <- as.phylo(treeNumbers, 8)
#' distances <- ClusteringInfoDist(trees)
#' clustering <- kmeans(distances, 4)
#' clust <- clustering$cluster
#' clust <- rep(1:4, c(21, 21, 21, 21))
#' nCl <- 4
#' cmds <- cmdscale(distances)
#' silmds <- cmds
#' 
#' layout4 <- rbind(cbind(rep(0, 4), rep(1, 4), 2:5),
#'              cbind(rep(6, 4), rep(7, 4), 8:11),
#'              cbind(rep(12, 4), rep(13, 4), 14:17))
#' layout3 <- rbind(cbind(rep(0, 3), rep(1, 3), 2:4),
#'              cbind(rep(5, 3), rep(6, 3), 7:9),
#'              cbind(rep(10, 3), rep(11, 3), 12:14))
#' layout(1 + if(nCl == 4) layout4 else layout3)
#' par(mar = rep(1.5, 4))
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
ClusterMDS <- function(clust, distances, proposal = cmdscale(distances),
                       maxIter = 1e6, stability = 42,
                       weight = c("sil" = 1,
                                  "ave" = 1,
                                  "max" = 1, # spacing between clusters
                                  "ks" = 1,
                                  "sor" = 1,
                                  "sov" = 1,
                                  "dfm" = 1,
                                  "mnn" = 1,
                                  "ste" = 1),
                       trace = FALSE) {
  distances <- as.dist(distances)
  treeID <- combn(attr(distances, "Size"), 2)
  nClust <- length(unique(clust))
  if (nClust != max(clust)) {
    stop("Cluster numbering must be consecutive")
  }
  Measuring <- function(x) {
    if (nClust < 2 && x %in% c("sil")) {
      FALSE
    } else {
      !is.na(weight[x]) && weight[x] != 0
    }
  }
  Score <- function(orig, new) sqrt(sum((orig - new) ^ 2))
  RatioScore <- function(orig, new) Score(orig / mean(orig), new / mean(new))
  scores <- 0 * weight
  
  if (Measuring("ks") || Measuring("ave")) {
    clustID <- matrix(clust[treeID], 2)
    clustPairs <- combn(nClust + 1, 2)
    clustPairs[2, ] <- clustPairs[2, ] - 1L
    clustPairKeys <- apply(clustPairs, 2, paste, collapse = "-")
    clustKey <- factor(apply(clustID, 2, 
                             function(x) paste(sort(x), collapse = "-")),
                       clustPairKeys)
    clustPair <- vapply(seq_along(clustPairKeys), function(pair) {
      clustKey == clustPairKeys[pair]
    }, logical(length(clustKey)))
    #normD <- distances / mean(distances)
    #clustDists <- apply(clustPair, 2, function(x) normD[x])
    clustDists <- apply(clustPair, 2, function(x) distances[x])
    ClusterStats <- function(d) {
      smry <- if(is.null(dim(d))) {
        vapply(d, summary, summary(1:3))
      } else {
        matrix(summary(d[, 1]), ncol = 1,
               dimnames = list(names(summary(1:3)), NULL))
      }

      rbind(
        med = smry["Median", ],
        mean = smry["Mean", ],
        # mad = vapply(d, mad, double(1)),
        # q1 = smry["1st Qu.", ],
        # q3 = smry["3rd Qu.", ],
        #iqr = smry["3rd Qu.", ] - smry["1st Qu.", ],
        min = smry["Min.", ] / smry["Median", ],
        max = smry["Max.", ] / smry["Median", ]
        #range = smry["Max.", ] - smry["Min.", ]
      )
    }
    clustZero <- ClusterStats(clustDists)
    colnames(clustZero) <- clustPairKeys
  }
  
  if (Measuring("sor") || Measuring("sov")) {
    hiDimMap <- suppressWarnings(
      cmdscale(distances, k = min(100L, attr(distances, "Size") - 1L))
    )
    sorZero <- SumOfRanges(hiDimMap, clust)
    sovZero <- SumOfVariances(hiDimMap, clust)
    dfmZero <- DistanceFromMedian(hiDimMap, clust)
  }
  
  if (Measuring("mnn")) {
    mnnZero <- MeanNN(distances, clust)
  }
  
  if (Measuring("ste")) {
    steZero <- MeanMSTEdge(distances, clust)
  }
  
  if (Measuring("sil")) {
    sil <- cluster::silhouette(clust, distances)
    silWidth <- sil[, "sil_width"]
  }
  bestScore <- Inf
  mover <- 0L
  accepted <- logical(maxIter)
  la <- integer(maxIter)
  mostStable <- 0
  log <- matrix(0, maxIter, length(weight),
                dimnames = list(NULL, names(scores)))
  lastAccepted <- c("any" = 1, "cl" = 1, "one" = 1, "ten" = 1, "sil" = 1)
  if (nClust < 2) {
    lastAccepted <- lastAccepted[setdiff(names(lastAccepted), "sil")]
  }
  strategy <- "any"
  
  
  range <- apply(proposal, 2, range)
  scale <- range[2, ] - range[1, ]
  
  cli_progress_bar("Rescaling", total = stability, .auto_close = FALSE)
  for (i in seq_len(maxIter)) {
    dprime <- dist(proposal)
    
    if (Measuring("ave") || Measuring("ks")) {
      #proposal.dist <- dprime / mean(dprime)
      #propDists <- apply(clustPair, 2, function(x) proposal.dist[x])
      propDists <- apply(clustPair, 2, function(x) dprime[x])
    }
    
    if (Measuring("ks")) {
      scores["ks"] <- median(unlist(suppressWarnings(mapply(
          ks.test, clustDists, propDists,
          MoreArgs = list(exact = FALSE, simulate.p.value = FALSE)
        )["statistic", ])))
    }
    
    if (Measuring("ave") || Measuring("max")) {
      propStats <- ClusterStats(propDists)
      
      if (Measuring("ave")) { 
        scores["ave"] <- RatioScore(clustZero[c("mean", "med"), ],
                                    propStats[c("mean", "med"), ])
      }
      if (Measuring("max")) {
        scores["max"] <- Score(clustZero[c("min", "max"), ],
                               propStats[c("min", "max"), ])
      }
    }
    
    if (Measuring("sor")) {
      scores["sor"] <- RatioScore(sorZero, SumOfRanges(proposal, clust))
    }
    
    if (Measuring("sov")) {
      scores["sov"] <- RatioScore(sovZero, SumOfVars(proposal, clust))
    }
    
    if (Measuring("dfm")) {
      scores["dfm"] <- RatioScore(dfmZero, DistFromMed(dprime, clust))
    }
    
    if (Measuring("mnn")) {
      scores["mnn"] <- RatioScore(mnnZero, MeanNN(dprime, clust))
    }
    
    if (Measuring("ste")) {
      scores["ste"] <- RatioScore(steZero, MeanMSTEdge(dprime, clust))
    }
    
    if (Measuring("sil")) {
      mappedSil <- cluster::silhouette(clust, dprime)
      dSil <- mappedSil[, "sil_width"] - silWidth
      dSil2 <- dSil * dSil
      d2Cum <- cumsum(dSil2)
      d2Sum <- d2Cum[length(d2Cum)]
      scores["sil"] <- d2Sum
    }
    
    score <- sum(scores * weight)
    
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
    
    strategyLottery <- cumsum(
      1 / lastAccepted[setdiff(names(lastAccepted), "any")]
    )
    switch(names(which.max(
      strategyLottery > runif(1, max = rev(strategyLottery)[1])
      )),
      "cl" = {
        # move cluster
        strategy <- "cl"
        mover <- clust == sample.int(nClust, 1)
      }, "ten" = {
        # move ten
        strategy <- "ten"
        mover <- sample.int(nrow(proposal), min(nrow(proposal), 10))
      }, "one" = {
        # move one
        strategy <- "one"
        mover <- sample.int(nrow(proposal), 1)
      }, "sil" = {
        # move one
        strategy <- "sil"
        mover <- which.max(d2Cum > runif(1, max = d2Sum))
      }
    )
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
