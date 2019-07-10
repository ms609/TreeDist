#' Mean distances between random pairs of trees
#' 
#' A three-dimensional array listing the normalized distances between 100 random 
#' pairs of trees drawn from the uniform distribution using
#' `ape::rtree(nTip, br=NULL)`.
#' 
#' Rows are named with an abbreviation of the tree comparison metric.
#' Columns list the mean and standard deviation of calculated tree distances.
#' The third dimension lists the number of tips in the trees compared.
#' 
#' @keywords datasets
"randomTreeDistances"

#' Distances between random pairs of 25-tip trees
#' 
#' A two-dimensional matrix listing the normalized distances between 10000 random
#' pairs of 25-tip trees drawn from the uniform distribution using
#' `ape::rtree(nTip, br=NULL)`.
#' 
#' Rows are named with an abbreviation of the tree comparison metric.
#' # TODO list these and state how normalized.
#' Columns list the calculated tree distances.
#' 
#' @keywords datasets
"distanceDistribution25"