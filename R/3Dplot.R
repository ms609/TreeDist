#' @export
#' @keywords internal
ThreeDPlotServer <- function(input, output, distances, clusterings,
                             mapping, mstEnds, mstSize,
                             pointCols, thinnedTrees,
                             StrainCol, mode3D) {
  
  if (!requireNamespace("rgl", quietly = TRUE)) {
    output$threeDPlot <- shiny::renderUI({
      shiny::div(
        class = "alert alert-warning",
        "3D visualisation requires the ",
        shiny::tags$code("rgl"),
        " package, which is not available in this environment."
      )
    })
    return(invisible(NULL))
  }
  
  rgl <- asNamespace("rgl")
  
  output$threeDPlot <- rgl$renderRglwidget({
    if (mode3D() && inherits(distances(), "dist")) {
      cl   <- clusterings()
      proj <- mapping()
      
      shiny::withProgress(message = "Drawing 3D plot", {
        rgl$open3d(useNULL = TRUE)
        shiny::incProgress(0.1)
        
        rgl$bg3d(color = "white")
        rgl$plot3d(
          proj[, 1], proj[, 2], proj[, 3],
          aspect = 1,
          axes   = FALSE,
          col    = pointCols(),
          alpha  = input$pt.opacity / 255,
          cex    = input$pt.cex,
          xlab = "", ylab = "", zlab = ""
        )
        
        shiny::incProgress(0.6)
        
        if ("labelTrees" %in% input$display) {
          rgl$text3d(proj[, 1], proj[, 2], proj[, 3], thinnedTrees())
        }
        
        if (mstSize() > 0) {
          rgl$segments3d(
            proj[t(mstEnds()), ],
            col = if ("mstStrain" %in% input$display) {
              rep(StrainCol(distances(), proj[, 1:3]), each = 2)
            } else {
              "#bbbbbb"
            },
            lty = 1
          )
        }
      })
      
      rgl$rglwidget()
    }
  })
}
