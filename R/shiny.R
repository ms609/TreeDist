#' Interactive tree space analysis
#' 
#' Explore tree landscapes in shiny app.
#' 
#' @template MRS
#' @family tree space functions
#' @export
TreeSpace <- function() {
  appDir <- system.file("treespace", package = "TreeDist")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing 'TreeDist'.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}