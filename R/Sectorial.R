#' Sectorial Search
#'
#' \code{SectorialSearch} performs a sectorial search on a tree, preserving the position of the root.
#'
#' \code{InapplicableSectorial} performs a sectorial search on the tree specified. A sectorial search 
#' detaches a random part of the tree, performs rearrangments on this subtree, then reattaches it 
#' to the main tree (Goloboff, 1999).
#' The improvement to local \var{pscore} hopefully (but not necessarily) improves the overall \var{pscore}.
#' As such, the output of \code{InapplicableSectorial} should be treated by further \acronym{TBR (/SPR/NNI)}
#' rearrangements and only retained if the ultimate parsimony score is better than 
#' that of the original tree.
#' 
#' \code{SectorialSearch} is a basic recipe that runs \code{InapplicableSectorial} followed by a few rounds
#' of tree rearrangement, returning a tree whose \var{pscore} is no worse than that of \code{start.tree}.
#' 
#' @param tree a rooted, resolved tree in \code{\link{phylo}} format from which to start the search;
#' @template datasetParam
#' @param outgroup a vector listing the taxa that form the outgroup;
#' @template concavityParam
#' @param maxit maximum number of sectorial iterations to perform;
#' @param maxIter maximum number of rearrangements to perform on each sectorial iteration;
#' @param cluster a cluster prepared using \code{\link{PrepareCluster}}; may speed up search on multicore machines;
#' @param k stop when \code{k} searches have improved their sectorial score;
#' @param verbosity integer determining how verbose the reporting to stdout will be;
#' @param smallest.sector sectors with fewer than \code{smallest.sector} taxa will not be selected; \kbd{4} is the smallest sensible value;
#' @param largest.sector sectors with more than \code{largest.sector} taxa will not be selected;
#' @param rearrangements method to use when rearranging subtrees: NNI, SPR or TBR;
#' @param \dots other arguments to pass to subsequent functions.
#' 
#' @return a rooted tree of class \code{phylo}.
#' 
#' @references Goloboff, P. (1999). \cite{Analyzing large data sets in reasonable times: solutions for composite optima.}
#' Cladistics, 15(4), 415-428. doi:\href{http://dx.doi.org/10.1006/clad.1999.0122}{10.1006/clad.1999.0122}
#' 
#' @author Martin R. Smith
#' 
#' @seealso \code{\link{TreeSearch}}
#' @seealso \code{\link{Ratchet}}
#' 
#' @examples{
#' data('SigSut')
#' outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
#' njtree <- Root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
#' njtree$edge.length <- NULL; njtree<-SetOutgroup(njtree, outgroup)
#' InapplicableSectorial(njtree, SigSut.phy, outgroup, maxit=1, maxIter=50, largest.sector=7)
#' \dontrun{SectorialSearch(njtree, SigSut.phy, outgroup, 'SPR') # Will be time-consuming}
#' 
#' ## SectorialSearch is currently defined as
#' # function (start.tree, dataset, outgroup, rearrangements='NNI') {
#' #   best.score <- attr(start.tree, 'pscore')
#' #   if (length(best.score) == 0) best.score <- InapplicableParsimony(start.tree, dataset)
#' #   sect <- InapplicableSectorial(start.tree, dataset, outgroup=outgroup, verbosity=0, maxit=30,
#' #               maxIter=200, maxHits=15, smallest.sector=6,
#  #               largest.sector=length(start.tree$edge[,2])*0.25, rearrangements=rearrangements)
#' #   sect <- DoTreeSearch(sect, dataset, outgroup, method='NNI', 
#' #                        maxIter=2000, maxHits=20, verbosity=3)
#' #   sect <- DoTreeSearch(sect, dataset, outgroup, method='TBR',
#' #                        maxIter=2000, maxHits=25, verbosity=3)
#' #   sect <- DoTreeSearch(sect, dataset, outgroup, method='SPR',
#' #                        maxIter=2000, maxHits=50, verbosity=3)
#' #   sect <- DoTreeSearch(sect, dataset, outgroup, method='NNI',
#' #                        maxIter=2000, maxHits=50, verbosity=3)
#' #   if (attr(sect, 'pscore') <= best.score) {
#' #     return (sect)
#' #   } else return (SetOutgroup(start.tree, outgroup))
#' # }
#' }
#' 
#' @keywords  tree 
#' @export
SectorialSearch <- function (tree, dataset, TreeScorer = FitchScore, concavity = NULL, rearrangements='NNI', maxIter=2000, cluster=NULL, verbosity=3, ...) {
  best.score <- attr(tree, 'pscore')
  tree <- RenumberTips(tree, names(dataset))
  if (length(best.score) == 0) best.score <- TreeScorer(tree, dataset, ...)[[1]]
  sect <- InapplicableSectorial(tree, dataset, cluster=cluster,
    verbosity=verbosity-1, maxit=30, maxIter=maxIter, maxHits=15, smallest.sector=6, 
    largest.sector=length(tree$edge[,2L])*0.25, rearrangements=rearrangements)
  sect <- DoTreeSearch(sect, dataset, method='NNI', maxIter=maxIter, maxHits=30, cluster=cluster, verbosity=verbosity)
  sect <- DoTreeSearch(sect, dataset, method='TBR', maxIter=maxIter, maxHits=20, cluster=cluster, verbosity=verbosity)
  sect <- DoTreeSearch(sect, dataset, method='SPR', maxIter=maxIter, maxHits=50, cluster=cluster, verbosity=verbosity)
  sect <- DoTreeSearch(sect, dataset, method='NNI', maxIter=maxIter, maxHits=60, cluster=cluster, verbosity=verbosity)
  if (attr(sect, 'pscore') <= best.score) {
    return (sect)
  } else return (tree)
}

#' Sectorial Search with inapplicable data
#'
#' \code{InapplicableSectorial} is called by SectorialSearch
#'
#' @param PARAM is a parameter you should send to it
#' 
#' @examples
#' to_do <- TRUE
#' 
#' @return This function returns :
#'   
#' @author Martin R. Smith
#' @importFrom ape root
#' @export
InapplicableSectorial <- function (tree, dataset, TreeScorer = FitchScore, maxit=100, 
    maxIter=500, k=5, verbosity=0, smallest.sector=4, largest.sector=1e+06, rearrangements="NNI", ...) {
  if (class(dataset) != 'phyDat') stop("dataset must be a phyDat object.")
  if (is.null(tree)) stop("a starting tree must be provided")
  tree <- RenumberTips(tree, names(dataset))
  if (verbosity >= 0) cat('InapplicableSectorial search: optimizing sectors of', smallest.sector, 'to', floor(largest.sector), 'tips')
  
  SectorData <- function (X, tips) {
    at <- attributes(X)
    dec <- X[tips, ]
    nBits <- floor(log2(max(X))) + 1L
    bin <- array(FALSE, dim=c(nrow(dec), ncol(dec), nBits))
    for (i in 0:(nBits-1)) {
      bin[, , nBits-i] <- as.logical(dec %% 2)
      dec <- (dec %/% 2)
    }
    state.union <- apply(bin, c(1,3), all)
    parsimony.informative <- !as.logical(rowSums(state.union))
    if (!any(parsimony.informative)) return (NULL)
    X <- X[tips, parsimony.informative]
    informative.chars <- sum(parsimony.informative)
    SECTOR_ROOT <- rep(2^nBits-1, informative.chars)
    X <- cbind(X, SECTOR_ROOT)
    attr(X, 'nr') <- informative.chars
    attr(X, 'inapp.level') <- at$inapp.level
    inapp.power2 <- log2(at$inapp.level) + 1
    #attr(X, 'min.steps') <- apply(X, 1, function(x) min.steps(x, inapp.power2))
    attr(X, 'levels') <- at$levels
    attr(X, 'weight') <- at$weight[parsimony.informative]
    class(X) <- 'morphyDat'
    warning("#TODO this is not yet tested")
    X
  }
  
  eps <- 1e-08
  kmax <- 1
  for (i in 1:maxit) {
    edge1 <- tree$edge[,1]
    nodes <- unique(edge1)[-1]
    node.lengths <- sapply(GetDescendants(tree, nodes), length) # (10x quicker than DoDescendants)
    candidate.nodes <- nodes[node.lengths >= smallest.sector & node.lengths <= largest.sector]
    if (verbosity >= 0) cat ("\n - Iteration", i, "- attempting sectorial search on node ")
    repeat {
      sector <- sample(candidate.nodes, 1)
      crown <- ExtractClade(tree, sector)
      crown.tips <- crown$tip.label
      sector.size <- length(crown.tips)
      cat(sector, 'size', sector.size, '...')
      crown.data <- SectorData(dataset, crown.tips)
      if (!is.null(crown.data)) break else cat('unsuitable (no dataset); trying')
      candidate.nodes <- candidate.nodes[-which(candidate.nodes==sector)]
      if (length(candidate.nodes == 0)) stop('No selectable sectors contain parsimony information! Either "largest.sector" is close to "smallest.sector" or your dataset is short of parsimony information.')
    } 
    if (verbosity >= 0) cat(' Sector OK.')
    crown <- root(AddTip(crown, 0, 'SECTOR_ROOT'), length(crown$tip.label) + 1, resolve.root=TRUE) ## TODO use Root or add ape::root to includeFrom in NAMESPACE
    initial.p <- TreeScorer(crown, crown.data, ...)
    attr(crown, 'pscore') <- initial.p
    if (verbosity >= 0) cat("\n - Running", rearrangements, "search on sector", sector)
    candidate <- TreeSearch(crown, crown.data, 'SECTOR_ROOT', method=rearrangements, verbosity=verbosity-1, maxIter=maxIter, ...)
    candidate.p <- attr(candidate, 'pscore')
    
    if((candidate.p + eps) < initial.p) {
      kmax <- kmax + 1
      stump <- DropTip(tree, GetDescendants(tree, sector)[[1]], subtree=TRUE)
      stump.edge <- 1:nrow(stump$edge)
      stump$root.edge <- 1
      crown <- DropTip(candidate, 'SECTOR_ROOT')
      tree <- CollapseSingles((BindTree(stump, crown, where=which(stump$tip.label==paste('[', sector.size, '_tips]', sep="")), position=0)))
      if (verbosity > 0) cat(' : improved local pscore, updated tree')
    } else if (verbosity > 0) cat (' : no improvement to local pscore')
    if (kmax == k) break()
  } # for
  if (verbosity >= 0)
    cat ("\nCompleted sectorial rearrangements.\n")
  attr(tree, 'pscore') <- NULL
  attr(tree, 'hits') <- NULL
  tree
}  # InapplicableSectorial

