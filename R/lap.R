#' Solve linear assignment problem using LAPJV
#'
#' Use the algorithm of \insertCite{Jonker1987;textual}{TreeDist} to solve the 
#' Linear Sum Assignment Problem (LSAP).
#' 
#' The Linear Assignment Problem seeks to match each row of a matrix with a 
#' column, such that the cost of the matching is minimized.
#' 
#' The Jonker & Volgenant approach is a faster alternative to the Hungarian
#' algorithm \insertCite{Munkres1957}{TreeDist}, which is implemented in 
#' `clue::solve_LSAP()`.
#' 
#' Note: the JV algorithm expects integers. In order to apply the function
#' to a non-integer _n_, as in the tree distance calculations in this package,
#' each _n_ is multiplied by the largest available integer before applying
#' the JV algorithm.  If two values of _n_ exhibit a trivial difference -- 
#' e.g. due to floating point errors -- then this can lead to interminable
#' run times.  (If numbers of the magnitude of billions differ only in their
#' last significant digit, then the JV algorithm may undergo billions of 
#' iterations.)  To avoid this, integers over 2^22 that differ by a value of
#' 8 or less are treated as equal.
#' 
#' @references \insertAllCited{}
#'
#' @author [C++ code](
#' https://github.com/yongyanghz/LAPJV-algorithm-c/blob/master/LAPJV/lap.cpp)
#' by Roy Jonker, MagicLogic Optimization Inc. <roy_jonker@magiclogic.com>, 
#' with contributions from Yong Yang <yongyanglink@gmail.com>, after 
#' [Yi Cao](https://uk.mathworks.com/matlabcentral/profile/authors/69713-yi-cao)
#' 
#' @param x Matrix of costs.
#' @return `LAPJV()` returns a list with two entries: `score`, the score of the
#' optimal matching;
#' and `matching`, the columns matched to each row of the matrix in turn.
#' 
#' @examples 
#' problem <- matrix(c(7, 9, 8, 9, 9,
#'                     2, 8, 5, 7, 9,
#'                     1, 6, 6, 9, 9,
#'                     3, 6, 2, 2, 9), 4, 5, byrow = TRUE)
#'
#' LAPJV(problem)
#' @seealso 
#' Implementations of the Hungarian algorithm exist in \pkg{adagio},
#' \pkg{RcppHungarian}, and \pkg{clue} and \pkg{lpSolve}; for larger matrices,
#' these are substantially slower. (See discussion at [Stack Overflow](
#' https://stackoverflow.com/questions/72806265/).)
#' 
#' The JV algorithm is implemented for square matrices in the Bioconductor
#' package [`GraphAlignment::LinearAssignment()`](
#' https://www.bioconductor.org/packages/release/bioc/html/GraphAlignment.html).
#' 
#' @export
LAPJV <- function(x) {
  dims <- dim(x)
  if (length(dims) == 2L) {
    if (any(dims == 0L)) {
      integer(0)
    } else {
      lapjv(x, max(x))
    }
  } else {
    stop("x must be a two-dimensional matrix")
  }
}
