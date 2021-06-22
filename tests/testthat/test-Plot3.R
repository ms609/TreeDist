test_that("Plot3() plots", {
  skip_if_not_installed('vdiffr', minimum_version = "1.0.0")
  library('vdiffr')
  skip_if(packageVersion("graphics") < "4.1")
  
  expect_doppelganger("Simple plot", {
    disorder <- c(1,5,2,6,3,7,4,8,5,9,10)
    Plot3(disorder, disorder, 6 - disorder, pch = 21, cex = 18,
          asp = 1, xlab = '', ylab = '',
          bg = 'white',
          plot.bg = '#bbeedd99',
          #shrink = 0,
          frame.plot = FALSE)
  })
  
  
  expect_doppelganger("Plotting order", {
    pts <- cbind(1:10, 11:20, c(1, 1, 1, 1, 10, 10, 10, 5, 6, 7))
    Plot3(pts,
          pch = c(3, rep(22, 4), rep(21, 4), 4),
          cex = c(rep(18, 6), 15, rep(18, 3)),
          col = c('red', rep(1, 4), rep(3, 4), 2),
          bg = c(rep('#ff555511', 5), rep('#bbbbff11', 5)),
          fog = 0.8, shrink = 0.8,
          asp = 1, xlab = '', ylab = '', frame.plot = FALSE)
  })
  
})