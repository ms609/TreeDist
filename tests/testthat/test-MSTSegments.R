test_that("MST example plots as expected", {
  skip_if_not_installed("graphics", "4.3")
  skip_if_not_installed("vdiffr", "1.0")
  vdiffr::expect_doppelganger("MST example plot", function() {
    distances <- structure(
      c(3, 2.3, 2.3, 2.3, 3, 1.7, 2.3, 2.3, 2.3, 4.5, 4.5, 1.7, 1.7, 2.3, 3,
        2.3, 1.7, 2.3, 2.3, 3.8, 4.4, 1.6, 3.1, 2.3, 2.6, 3, 1.5, 2.6, 4.7,
        4.6, 2.6, 2.3, 1.5, 3, 2.6, 1.5, 4.3, 4.7, 1.7, 2.6, 1.5, 3, 1.6, 4.6,
        4.7, 2.3, 2.3, 1.7, 1.7, 4.4, 3.8, 3.1, 3.1, 1.5, 4.4, 4.4, 3.1, 2.6,
        4.2, 4.8, 3, 4.8, 4.2, 4.7, 4.3, 1.6),
      class = "dist", Size = 12L, Diag = FALSE, Upper = FALSE
    ) # dput(round(ClusteringInfoDist(as.phylo(5:16, 8)), 1))
    
    mapping <- structure(
      c(0.75, 0.36, 0.89, 0.85, 0.89, 0.36, 0.64, 0.6, 0.6, 0.85, -3.4, -3.4,
        0, -0.25, 1.33, 0.39, -1.33, 0.25, 0, -1.52, 1.52, -0.39, -0.58, 0.58),
      dim = c(12L, 2L)
    ) # dput(round(cmdscale(distances, k = 2), 2))
    
    mstEnds <- MSTEdges(distances)
    
    # Set up blank plot
    plot(mapping, asp = 1, frame.plot = FALSE, ann = FALSE, axes = FALSE,
         type = "n")
    
    # Add MST
    MSTSegments(mapping, mstEnds, 
                col = StrainCol(distances, mapping, mstEnds))
    
    # Add points at end so they overprint the MST
    points(mapping)
    PlotTools::SpectrumLegend(
      "bottomleft",
      legend = c("Extended", "Median", "Contracted"),
      bty = "n",
      y.intersp = 2,
      lend = "square",
      palette = hcl.colors(256L, "RdYlBu", rev = TRUE)
    )
  })
})

test_that("StrainCol() handles zeroes", {
  distances <- dist(c(1, 1, 10, 100))
  mapping <- cmdscale(distances)
  expect_equal(attr(StrainCol(distances, mapping), "logStrain")[1], Inf)
})
