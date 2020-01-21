library('TreeTools')

#' @noRd
OneOverlap <- function(A1, A2, nLeaves) {
  if (A1 == A2) {
    # Return:
    LnRooted.int(A1) + LnRooted.int(nLeaves - A2)
  } else {
    if (A1 < A2) {
      # Return:
      LnRooted.int(A2) + LnRooted.int(nLeaves - A1) - LnRooted.int(A2 - A1 + 1L) 
    } else {
      # Return:
      LnRooted.int(A1) + LnRooted.int(nLeaves - A2) - LnRooted.int(A1 - A2 + 1L) 
    }
  }
}

oneOverlap <- lapply(seq_len(100), function (n) 
  matrix(mapply(OneOverlap, rep(seq_len(n), n), rep(seq_len(n), each=n),
                nLeaves = n), ncol=n))


usethis::use_data(oneOverlap, compress='xz', overwrite=TRUE, internal=TRUE)
