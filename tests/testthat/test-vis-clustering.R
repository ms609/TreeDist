test_that("ClusterMDS()", {
  
  library("TreeTools", quietly = TRUE)
  treeNumbers <- c(0:10, 40:50, 120:130, 200:210)
  trees <- as.phylo(treeNumbers, 8)
  distances <- ClusteringInfoDist(trees)
  clustering <- kmeans(distances, 4)
  clust <- clustering$cluster
  cmds <- cmdscale(distances)
  
  par(mfrow = c(2, 3))
  plot(cluster::silhouette(clust, distances), col = 1:4, main = "Original")
  
  silmds <- ClusterMDS(clust, distances, cmds, stability = 42)
  plot(silmds, col = clust, cex = 1.5, asp = 1, axes = FALSE, ann = FALSE)
  expect_equal(cluster::silhouette(clust, dist(silmds))[, "sil_width"],
               cluster::silhouette(clust, distances)[, "sil_width"],
               tolerance = 0.01)
  
})
