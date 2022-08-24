points <- rbind(matrix(1:16, 4),
                rep(1, 4),
                matrix(1:32, 8, 4) / 10)
cluster <- rep(1:3, c(4, 1, 8))

test_that("SumOfRanges()", {
  expect_warning(expect_equal(SumOfRanges(points[cluster == 2, ]), 0),
                 "lacks dimensions")
  expect_equal(SumOfRanges(points, cluster), c(3 * 4, 0, 0.7 * 4))
  expect_equal(
    sapply(1:3, function(i) SumOfRanges(points[cluster == i, , drop = FALSE])),
    SumOfRanges(points, cluster)
  )
})

test_that("SumOfVariances()", {
  expect_warning(expect_equal(SumOfVariances(points[cluster == 2, ]), NA_real_),
                 "lacks dimensions")
  expect_equal(SumOfVars(points, cluster),
               c(var(1:4) * 4, NA_real_, var(1:8 / 10) * 4))
  expect_equal(
    sapply(1:3, function(i) SumOfVars(points[cluster == i, , drop = FALSE])),
    SumOfVariances(points, cluster)
  )
})

test_that("MeanCentroidDistance()", {
  expect_warning(expect_equal(MeanCentroidDistance(points[cluster == 2, ]), 0),
                 "lacks dimensions")
  
  expect_equal(MeanCentroidDistance(points, cluster), c(2, 0, 0.4))
  expect_equal(
    sapply(1:3, function(i) MeanCentDist(points[cluster == i, , drop = FALSE])),
    MeanCentroidDist(points, cluster)
  )
})

test_that("MeanNN()", {
  expect_warning(expect_equal(MeanNN(points[cluster == 2, ]), NA_real_),
                 "lacks dimensions")
  expect_equal(MeanNN(points, cluster), c(2, NA, 0.2))
  expect_equal(
    sapply(1:3, function(i) MeanNN(points[cluster == i, , drop = FALSE])),
    MeanNN(points, cluster)
  )
  expect_equal(MeanNN(points, cluster), MeanNN(dist(points), cluster))
})

test_that("MeanMSTEdge()", {
  expect_warning(expect_equal(MeanMSTEdge(points[cluster == 2, ]), NA_real_),
                 "lacks dimensions")
  expect_equal(MeanMSTEdge(points, cluster),
               c(2, NA, 0.2))
  expect_equal(
    sapply(1:3, function(i) MeanMSTEdge(points[cluster == i, , drop = FALSE])),
    MeanMSTEdge(points, cluster)
  )
  expect_equal(MeanMSTEdge(dist(points), cluster), MeanMSTEdge(points, cluster))
})
