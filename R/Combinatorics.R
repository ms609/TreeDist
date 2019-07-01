globalVariables(c('doubleFactorials', 'logDoubleFactorials'), 'TreeSearch')

#' Double Factorial
#' 
#' @param n Vector of integers.
#' 
#' @return Returns the double factorial, n x (n - 2) x (n - 4) x (n - 6) x ...
#' 
#' @examples {
#' DoubleFactorial (-4:0) # Return 1 if n < 2
#' DoubleFactorial (2) # 2
#' DoubleFactorial (5) # 1 x 3 x 5
#' exp(LogDoubleFactorial.int (8)) # 2 x 4 x 6 x 8
#' 
#' }
#' 
#' @author Martin R. Smith
#' @family Double factorial
#' @export
DoubleFactorial <- function (n) {
  if (any(n > 300)) stop("301!! is too large to store as an integer. Use LogDoubleFactorial instead.")
  
  n[n < 2] <- 1
  doubleFactorials[n]
  
  #
  #odds <- as.logical(x %% 2)
  #
  #oddX <- x[odds]
  #xPlusOneOverTwo <- (oddX + 1) / 2
  #evenX <- x[!odds]
  #xOverTwo <- evenX / 2
  #
  #ret <- integer(length(x))
  #ret[odds] <- gamma(oddX + 1L) / (gamma(xPlusOneOverTwo) * 2^(xPlusOneOverTwo - 1L))
  #ret[!odds] <- evenX * gamma(xOverTwo) * 2^(xOverTwo - 1L)
  #
  ## Return:
  #ret
}

# Memoizing this function makes it MUCH slower...
#' @describeIn DoubleFactorial Returns the logarithm of the double factorial.
LogDoubleFactorial <- (function (n) {
  n[n < 2] <- 1 # Much faster than pmax
  if (any(n > 49999L)) {
    
    odds <- as.logical(n %% 2)
    
    oddN <- n[odds]
    nPlusOneOverTwo <- (oddN + 1) / 2
    evenN <- n[!odds]
    nOverTwo <- evenN / 2
    
    ret <- integer(length(n))
    ret[odds] <- lgamma(oddN + 1L) - (lgamma(nPlusOneOverTwo) + (nPlusOneOverTwo - 1) * log(2))
    ret[!odds] <- log(evenN) + lgamma(nOverTwo) + (nOverTwo - 1) * log(2)
    
    # Return:
    ret
    
  } else {
    # Return from cache
    logDoubleFactorials[n]
  }
})

#' @describeIn DoubleFactorial Slightly faster, when x is known to be length one
#' and below 50001
#' @export
LogDoubleFactorial.int <- function (n) {
  if (n < 2) {
    0
  } else {
    logDoubleFactorials[n]
  }
}

#' Number of rooted/unrooted trees
#' These functions return the number of rooted or unrooted trees consistent with a given pattern
#'  of splits.
#' 
#' Functions starting N return the number of rooted or unrooted trees, functions starting Ln
#' provide the log of this number.  Calculations follow Carter et al. 1990, Theorem 2.
#'
#' @param tips The number of tips.
#' @param splits vector listing the number of taxa in each tree bipartition.
#'
#' @author Martin R. Smith
#' 
#' @references 
#'  \insertRef{Carter1990}{TreeSearch}
#'  
#' @examples
#'   NRooted(10)
#'   NUnrooted(10)
#'   LnRooted(10)
#'   LnUnrooted(10)
#'   # Number of trees consistent with a character whose states are 00000 11111 222
#'   NUnrootedMult(c(5,5,3))
#' 
#' @export
NRooted     <- function (tips)  DoubleFactorial(tips + tips - 3L) # addition faster than 2*
#' @describeIn NRooted Number of unrooted trees
#' @export
NUnrooted   <- function (tips)  DoubleFactorial(tips + tips - 5L)
#' @describeIn NRooted  Log Number of unrooted trees
#' @export
LnUnrooted  <- function (tips) LogDoubleFactorial(tips + tips - 5L)
#' @describeIn NRooted  Log Number of unrooted trees (as integer)
#' @export
LnUnrooted.int <- function (tips) if (tips < 3) 0 else {
  logDoubleFactorials[tips + tips - 5L]
}
#' @describeIn NRooted  Log Number of rooted trees
#' @export
LnRooted    <- function (tips) LogDoubleFactorial(tips + tips - 3L)
#' @describeIn NRooted  Log Number of rooted trees (as integer)
#' @export
LnRooted.int <- function (tips) {
  if (tips < 2L) 0 else logDoubleFactorials[tips + tips - 3L]
}

#' Number of trees one SPR step away
#' Formula given by Allen and Steel (2001).
#'
#' @param n Number of tips in tree.
#' @references 
#'  \insertRef{Allen2001}{TreeSearch}
#' 
#' @export
N1Spr <- function (n) if (n > 3L) (n + n - 6L) * (n + n - 7L) else 0L
