library(TreeSearch)

#' @noRd
OneOverlap <- function(A1, A2, nTerminals) {
  if (A1 == A2) {
    # Return:
    LnRooted.int(A1) + LnRooted.int(nTerminals - A2)
  } else {
    if (A1 < A2) {
      # Return:
      LnRooted.int(A2) + LnRooted.int(nTerminals - A1) - LnRooted.int(A2 - A1 + 1L) 
    } else {
      # Return:
      LnRooted.int(A1) + LnRooted.int(nTerminals - A2) - LnRooted.int(A1 - A2 + 1L) 
    }
  }
}

oneOverlap <- lapply(seq_len(100), function (n) 
  matrix(mapply(OneOverlap, rep(seq_len(n), n), rep(seq_len(n), each=n),
                nTerminals = n), ncol=n))


usethis::use_data(oneOverlap, compress='xz', overwrite=TRUE, internal=TRUE)
