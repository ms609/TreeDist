#' Solve linear assignment problem using LAPJV
#'
#' Solve the 
#' [Linear Sum Assignment Problem](http://www.assignmentproblems.com/doc/LSAPIntroduction.pdf)
#' with the algorithm of Jonker & Volgenant (1987).
#' 
#' The Linear Assignment Problem seeks to match each row of a matrix with a 
#' column, such that the cost of the matching is minimized.
#' 
#' The JV approach improves on the Hungarian algorithm, which is implemented 
#' in [`clue::solve_LSAP`].
#' 
#' 
#' @references 
#' 
#' \insertRef{Jonker1987}{TreeDist}
#'
#' @author [C++ code](https://github.com/yongyanghz/LAPJV-algorithm-c/blob/master/LAPJV/lap.cpp)
#' by Roy Jonker, MagicLogic Optimization Inc. <roy_jonker@magiclogic.com>, 
#' with contributions from Yong Yang <yongyanglink@gmail.com>, after 
#' [Yi Cao](https://uk.mathworks.com/matlabcentral/profile/authors/69713-yi-cao)
#' 
#' 
#' @param x Square matrix of costs.
#' @return A list with two entries: `score`, the score of the optimal matching;
#' and `matching`, the columns matched to each row of the matrix in turn.
#' 
#' @examples 
#' 
#' problem <- matrix(c(7, 9, 8, 9,
#'                     2, 8, 5, 7,
#'                     1, 6, 6, 9,
#'                     3, 6, 2, 2), 4, 4, byrow=TRUE)
#'
#' LAPJV(problem)
#' @export
LAPJV <- function (x) {
  lapjv(x, max(x))
}