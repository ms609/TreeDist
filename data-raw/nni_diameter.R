library('TreeTools')

TreeNumbers <- function (nTip) {
  shapes <- lapply(seq_len(NUnrootedShapes(nTip)) - 1L, UnrootedTreeWithShape, nTip)
  numbers <- vapply(shapes, as.TreeNumber, double(1))
  paste0(".\\NNI-diameter-", nTip, ".exe ", numbers, collapse="; ")
}

TreeNumbers(12L)
