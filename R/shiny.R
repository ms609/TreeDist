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
  if (!requireNamespace('cluster', quietly = TRUE)) install.packages('cluster')
  if (!requireNamespace('protoclust', quietly = TRUE)) install.packages('protoclust')
  if (!requireNamespace('MASS', quietly = TRUE)) install.packages('MASS')
  if (!requireNamespace('Quartet', quietly = TRUE)) install.packages('Quartet')
  if (!requireNamespace('rgl', quietly = TRUE)) install.packages('rgl')
  if (!requireNamespace('readxl', quietly = TRUE)) install.packages('readxl')
  if (!requireNamespace('viridisLite', quietly = TRUE)) install.packages('viridisLite')
  
  shiny::runApp(appDir, display.mode = "normal")
}
