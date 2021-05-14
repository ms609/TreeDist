test_that("Plot3() plots", {
  skip_if_not_installed('vdiffr')
  library('vdiffr')
  skip_if(packageVersion("graphics") > "4.0.99")
  
  expect_doppelganger("Simple plot", {
    disorder <- c(1,5,2,6,3,7,4,8,5,9,10)
    Plot3(disorder, disorder, 6 - disorder, pch = 21, cex = 15,
          asp = 1, xlab = '', ylab = '',
          bg = 'white',
          plot.bg = '#bbeedd99',
          
          frame.plot = FALSE)
  })
  
  
  expect_doppelganger("Plotting order", {
    Plot3(1:10, 11:20, c(1, 1, 1, 1, 10, 10, 10, 5, 6, 7),
          pch = c(3, rep(22, 4), rep(21, 4), 4),
          cex = c(rep(10, 6), 8, rep(10, 3)),
          col = c('red', rep(1, 4), rep(3, 4), 2),
          bg = c(rep('#ff555511', 5), rep('#bbbbff11', 5)),
          fog = 0.9, shrink = 0.9,
          asp = 1, xlab = '', ylab = '', frame.plot = FALSE)
  })
  
})