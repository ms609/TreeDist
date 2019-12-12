#' Solve linear assignment problem using LAPJV
#'
#' Solve the Linear Assignment Problem with the algorithm of 
#' Jonker & Volgenant (1987).
#' 
#' The Linear Assignment Problem seeks to match each row of a matrix with a 
#' column, such that the cost of the matching is minimized.
#' 
#' The JV approach improves on the Hungarian algorithm, which is implemented 
#' in [`clue::solve_LSAP`].
#' 
#' 
#' @references 
#' "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear
#'  Assignment Problems," Computing 38, 325-340, 1987
#' 
#' R. Jonker and A. Volgenant, University of Amsterdam.
#' 
#'
#' @author C++ code by Roy Jonker @ MagicLogic Optimization Inc. 
#' e-mail: roy_jonker@magiclogic.com
#' 
#' Changed 2016-05-13 by Yong Yang(yongyanglink@gmail.com) in column reduction 
#' part according to matlab version of LAPJV algorithm (Yi Cao)
#' https://github.com/yongyanghz/LAPJV-algorithm-c/blob/master/LAPJV/lap.cpp
#' 
#' @param x Matrix to solve.
#' @return A list with two entries: `score`, the score of the optimal matching;
#' and `matching`, the columns matched to each row of the matrix in turn.
#' @export
LAPJV <- function (x) {
  lapjv(x, max(x))
}