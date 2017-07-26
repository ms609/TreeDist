#' @importFrom ape consensus
#' @export
SuccessiveApproximations <- function (tree, data, outgroup = NULL, k = 3, max.succiter = 20,
                                      pratchhits = 100, searchhits = 50, searchiter = 500,
                                      pratchiter = 5000, track = 0, suboptimal = 0.1) {
  
  if (k < 1) stop ('k should be at least 1, see Farris 1969 p.379')
  attr(data, 'sa.weights') <- rep(1, length(attr(data, 'weight')))
  collect.suboptimal <- suboptimal > 0
  
  max.node <- max(tree$edge[, 1])
  n.tip <- length(tree$tip.label)
  n.node <- max.node - n.tip
  bests <- vector('list', max.succiter + 1)
  bests.consensus <- vector('list', max.succiter + 1)
  best <- bests[[1]] <- bests.consensus[[1]] <- Root(tree, outgroup) #TODO does ape::root work here?
  for (i in seq_len(max.succiter) + 1) {
    if (track > 0) cat('\nSuccessive Approximations Iteration', i - 1)
    attr(best, 'score') <- NULL
    if (suboptimal > 0) {
      suboptimal.search <- suboptimal * sum(attr(data, 'sa.weights') * attr(data, 'weight'))
    }
    trees <- Ratchet(best, data, ParsimonyScorer = SuccessiveWeights, all = collect.suboptimal, 
                           suboptimal=suboptimal.search,    rearrangements='NNI',
                           pratchhits=pratchhits, searchhits=searchhits, searchiter=searchiter, 
                           pratchiter=pratchiter, outgroup = outgroup, track=track - 1)
    trees <- unique(trees)
    bests[[i]] <- trees
    suboptimality <- Suboptimality(trees)
    bests.consensus[[i]] <- ape::consensus(trees[suboptimality == 0])
    if (all.equal(bests.consensus[[i]], bests.consensus[[i - 1]])) return(bests[2:i])
    best <- trees[suboptimality == 0][[1]]
    l.i <- Fitch(best, data)
    p.i <- l.i / (n.node - 1)
    w.i <- ((p.i)^-k) - 1
    attr(data, 'sa.weights') <- w.i
  }
  cat('Stability not reached.')
  return(bests)
}

#' @keywords internal
#' @export
Suboptimality <- function (trees, proportional = FALSE) {
  scores <- vapply(trees, attr, double(1), 'score')
  if (proportional) {
    return ((scores - min(scores)) / min(scores))
  } else {
    return(scores - min(scores))
  }
}

#' @keywords internal
#' @export
SuccessiveWeights <- function(tree, data) {
  # Data
  if (class(data) == 'phyDat') data <- PrepareDataSA(data)
  if (class(data) != 'saDat') {
    stop('Invalid data type; prepare data with PhyDat() or PrepareDataSA().')
  }
  at <- attributes(data)
  weight <- at$weight
  sa.weights <- at$sa.weights
  if (is.null(sa.weights)) sa.weights <- rep(1, length(weight))
  steps <- Fitch(tree, data, at)
  return(sum(steps * sa.weights * weight))
}

PrepareDataSA <- function (data) {
# Written with reference to phangorn:::prepareDataFitch
  at <- attributes(data)
  nam <- at$names
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(data, "contrast")
  nTip <- length(data)
  
  at$names <- NULL
  powers.of.2 <- 2L ^ c(0L:(nLevel - 1L))
  tmp <- cont %*% powers.of.2
  tmp <- as.integer(tmp)
  data <- unlist(data, FALSE, FALSE)
  ret <- tmp[data] 
  ret <- as.integer(ret)
  attributes(ret) <- at
  inappLevel <- which(at$levels == "-")
  attr(ret, 'inappLevel') <- 2 ^ (inappLevel - 1)
  attr(ret, 'dim') <- c(nChar, nTip)  
  attr(ret, 'unique.tokens') <- apply(ret, 1, function(x) min(x, inappLevel))
  applicableTokens <- setdiff(powers.of.2, 2 ^ (inappLevel - 1))
  attr(ret, 'split.sizes') <- t(apply(ret, 1, function(x) vapply(applicableTokens, function (y) sum(x == y), integer(1))))
  dimnames(ret) <- list(NULL, nam)
  class(ret) <- 'saDat'
  ret
}