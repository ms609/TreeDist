# Memoizing this function makes it MUCH slower...
#' @describeIn DoubleFactorial Returns the logarithm of the double factorial.
LogDoubleFactorial <- (function (x) {
  x[x < 2] <- 1
  odds <- as.logical(x %% 2)
  
  oddX <- x[odds]
  xPlusOneOverTwo <- (oddX + 1) / 2
  evenX <- x[!odds]
  xOverTwo <- evenX / 2
  
  ret <- integer(length(x))
  ret[odds] <- lgamma(oddX + 1L) - (lgamma(xPlusOneOverTwo) + (xPlusOneOverTwo - 1) * log(2))
  ret[!odds] <- log(evenX) + lgamma(xOverTwo) + (xOverTwo - 1) * log(2)
  
  # Return:
  ret
})

logDoubleFactorials <- vapply(seq_len(50000), LogDoubleFactorial, double(1))

DoubleFactorial <- function (x) {
  if (any(x > 300)) stop("301!! is too large to store as an integer. Use LogDoubleFactorial instead.")
  
  
  x[x < 2] <- 1
  odds <- as.logical(x %% 2)
  
  oddX <- x[odds]
  xPlusOneOverTwo <- (oddX + 1) / 2
  evenX <- x[!odds]
  xOverTwo <- evenX / 2
  
  ret <- integer(length(x))
  ret[odds] <- gamma(oddX + 1L) / (gamma(xPlusOneOverTwo) * 2^(xPlusOneOverTwo - 1L))
  ret[!odds] <- evenX * gamma(xOverTwo) * 2^(xOverTwo - 1L)
  
  # Return:
  ret
}

doubleFactorials <- exp(logDoubleFactorials[seq_len(300)]) # Greater than 300 -> "Inf"
usethis::use_data(logDoubleFactorials, doubleFactorials, overwrite=TRUE)
