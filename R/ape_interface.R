#' Ape Time
#' 
#' Reads the time that an ape tree was modified from the comment in the nexus file
#'
#' @param filename Character string specifying path to the file
#' @param format Format in which to return the time: 'double' as a sortable numeric; 
#'               any other value to return a string in the format YYYY-MM-DD hh:mm:ss
#'
#' @return The time that the specified file was created by ape.
#' @export
#' @author Martin R. Smith
#'
ApeTime <- function (filename, format='double') {
  comment <- readLines(filename, n=2)[2]
  Month <- function (month) {
    months <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
    whichMonth <- months == month
    if (any(whichMonth)) {
      formatC(which (whichMonth), width=2, flag="0")
    } else {
      month
    }
  }
  DATEEXP <- ".*? (\\w+)\\s(\\d+)\\s(\\d+\\:\\d\\d\\:\\d\\d)\\s(\\d\\d\\d\\d).*"
  time <- paste0(gsub(DATEEXP, "\\4-", comment),
                 Month(gsub(DATEEXP, "\\1", comment)),
                 gsub(DATEEXP, "-\\2 \\3", comment))
  
  # Return:
  ifelse(format=='double', as.numeric(as.POSIXct(time, tz = "GMT")), time)
}
