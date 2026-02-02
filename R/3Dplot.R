#' @export
#' @keywords internal
ThreeDPlotServer <- function(input, output, distances, clusterings,
                             mapping, mstEnds, mstSize,
                             pointCols, thinnedTrees,
                             StrainCol, mode3D) {
  output$threeDPlot <- shiny::renderUI({
    plotSize <- sprintf( "height: %spx; width: %spx;", 600, 600)
    
    if (!mode3D()) {
      return(shiny::div(style = plotSize))
    }
    
    if (!requireNamespace("rgl", quietly = TRUE)) {
      is_wasm <- identical(Sys.getenv("R_WASM"), "1")
      if (!is_wasm) {
        if (!requireNamespace("rgl", quietly = TRUE)) {
          tryCatch(install.packages("rgl"), error = function(e) {})
        }
      }
    }
    
    if (!requireNamespace("rgl", quietly = TRUE)) {
      msg <- if (is_wasm) {
        shiny::p(
          shiny::strong("3D plot unavailable."),
          shiny::br(),
          "3D plots are not supported in the browser version."
        )
      } else {
        shiny::p(
          shiny::strong("3D plot unavailable."),
          shiny::br(),
          "The ",
          shiny::tags$code("rgl"),
          " package is not installed."
        )
      }
      
      return(shiny::div(
        style = paste(plotSize, "
            display: flex;
            align-items: center;
            justify-content: center;
            border: 1px solid #ddd;
            background: #fafafa;
            color: #666;
          "),
        class = "alert alert-warning",
        msg
      ))
    }
    
    rgl <- asNamespace("rgl")
    
    return(rgl$renderRglwidget({
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
    }))
  })
}
