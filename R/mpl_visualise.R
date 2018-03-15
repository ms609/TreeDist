### #' Visualize Inheritance
### #'
### #' \code{VisualiseInheritance} does something useful
### #'
### #' @param PARAM is a parameter you should send to it
### #' 
### #' @examples
### #' to_do <- TRUE
### #' 
### #' @return This function returns :
### #'   
### #' @author Martin R. Smith
### #'
### #' @importFrom graphics par plot
### #' @export
### VisualiseInheritance <- VisualizeInheritance <- VisIn <- function (tree, data, char.no, plot.fun=plot) {
###   par(mfrow=c(1,2), mar=rep(0.5,4))
###   VisualizeCharacter(tree, data, char.no, plot.fun)
### }
### #' VisualiseCharacter
### #' 
### #' Visualize optimization for a character
### #' 
### #' @description Determine and depict the possible states for a character on a tree under the most parsimonious conditions
### #' @usage VisualiseCharacter(tree, dataset, char.no, plot.fun = plot)
### #' 
### #' @param tree a fully-resolved tree in \code{\link{phylo}} format, with the desired outgroup; edge lengths are not supported and will be deleted;
### #' @template datasetParam
### #' @param char.no number of the character to be displayed;
### #' @param plot.fun a function that plots a tree, \code{\link{plot}} by default.
### #'
### #' @return The function plots a cladogram, annotating each node with the values that could occur there at the lowest parsimony score.
### #'   Nodes highlighed in red potentially (but do not necessarily) contribute to tree length.
### #' 
### #' @author Martin R. Smith
### #' 
### #' @examples{
### #' data('SigSut')
### #' outgroup <- c('Lingula', 'Mickwitzia', 'Neocrania')
### #' njtree <- root(nj(dist.hamming(SigSut.phy)), outgroup, resolve.root=TRUE)
### #' njtree$edge.length <- NULL; njtree<-SetOutgroup(njtree, outgroup)
### #' VisualiseCharacter(njtree, SigSut.phy, 5, function(x) {plot(x); edgelabels();})
### #' }
### #' @keywords  tree 
### #' @aliases VisualiseCharacter
### #' @aliases VisualiseChar
### #' @aliases VisualizeCharacter
### #' @aliases VisualizeChar
### #' 
### #' @importFrom ape .PlotPhyloEnv tiplabels nodelabels 
### #' @importFrom graphics text
### #' @export
### VisualizeCharacter <- VisualiseCharacter <- VisualiseChar <- VisualizeChar <- 
### function (tree, dataset, char.no, plot.fun = plot) {
###   if (class(dataset) == 'phyDat') dataset <- MorphyDat(dataset)
###   if (class(dataset) != 'morphyDat') stop('Invalid dataset type in VizualizeCharacter.')
###   warning("#TODO: Update to use new morphyDat dataset objects")
###   at <- attributes(dataset)
###   if (char.no > at$nr || char.no < 1) stop(paste0("char.no must be between 1 and ", at$nr, ' (', sum(at$weight), 'TS, ', at$nr, ' unique)'))
###   char.dat <- dataset[char.no,]
###   char.index <- at$index[char.no]
###   if (is.null(at$order) || at$order != "postorder") tree <- Postorder(tree)
###   tree.edge <- tree$edge
###   parent <- tree.edge[,1]
###   child <- tree.edge[,2]
###   tip.label <- tree$tip.label
###   nEdge <- length(parent)
###   nTip <- length(tip.label)
###   nNode <- nTip - 1
###   maxNode <- nNode + nTip
###   inapp <- at$inapp.level
###   tips <- seq(nTip)
###   nodes <- nTip + seq(nNode)
###   parentof <- parent[match((nTip + 2L):maxNode, child )] # Exclude the root, which has no parent
###   childof <- child [c(match(nodes, parent), length(parent) + 1L - match(nodes, rev(parent)))]
###   if (any(is.na(char.dat[tip.label]))) stop("Tree's tip labels could not all be found in dataset matrix")
###   
###   plot.fun(tree)
###   ret <- .Call("FITCHINAPP", char.dat[tip.label], as.integer(1), as.integer(parent), as.integer(child), as.integer(parentof), as.integer(childof), as.integer(nEdge), as.integer(nNode), as.double(1), as.integer(maxNode), as.integer(nTip), as.integer(inapp), PACKAGE='inapplicable')
###   downpass.states <- ret[[3]]
###   down.scorers <- ret[[4]]
###   inapp.nodes <- ret[[5]] > 0
###   
###   down.change <- sapply(nodes, function(n) {
###       children <- child[parent==n]
###       return (down.scorers[n] != down.scorers[children[1]] + down.scorers[children[2]])   
###     })
###   text(1,1,paste0('Char ', char.no, ' - TS', paste(which(at$index == char.no), collapse=', '), ': +', ret[[1]]), pos=4, cex=0.8)
###   
###   
###   tipcols = c('#fafafa', '#fafafa', '#fafabb', '#ffbbbb', '#bbffbb', '#bbbbff', '#bbbbff', '#bbffbb', '#ffbbbb', '#bbddff', '#ffbbdd')
###   names(tipcols) <- c(NA, max(downpass.states[1,]), max(downpass.states[1,])-inapp, 2^(0:7))
###   tipcols[as.character(inapp)] <- '#999999'
###   tipcols <- rev(tipcols)
###   bgcols <- tipcols[as.character(downpass.states[1,tips])]
###   bgcols[is.na(bgcols)] <- '#ffffbb'
###   tiplabels(PossibleTokens(at$levels, downpass.states[1,tips]), adj=c(0.3,0.5), bg=bgcols, col='#000088', cex=0.85)
###   nodelabels(PossibleTokens(at$levels, downpass.states[1,nodes]), adj=rep(1.25,2), bg=tipcols[as.character(downpass.states[1,nodes])], font=ifelse(down.change, 2, 1) , col=ifelse(down.change, '#cc3333', '#880000cc'), cex=ifelse(down.change,1,0.6))
###   
###   nodelabels(ifelse(inapp.nodes[nodes], '+', '-'), adj=c(1.25,-0.75), col=ifelse(inapp.nodes[nodes], '#008800', '#880000'), frame='none')
### }
### #' TITLE GOES HERE
### #'
### #' \code{FUNCTIONNAME} does something useful
### #'
### #' @param PARAM is a parameter you should send to it
### #' 
### #' @examples
### #' to_do <- TRUE
### #' 
### #' @return This function returns :
### #'   
### #' @author Martin R. Smith
### #' @export
### PossibleTokens <- function (lvls, number) {
###   nTokens <- length(lvls)
###   nNumber <- length(number)
###   output <- function (x) {paste0(x, collapse='')}
###   if (nNumber == 1) {    
###     if (number == 2^nTokens - 1) return('?')
###     which.levels <- rep(FALSE, nTokens)
###     binary <- AsBinary(number)
###     which.levels[seq_along(binary)] <- binary
###     return (output(lvls[as.logical(which.levels)]))
###   }
###   which.levels <- matrix(FALSE, nNumber, nTokens)
###   binary <- AsBinary(number)
###   which.levels[,seq_along(binary[1,])] <- as.logical(binary)
###   apply(which.levels, 1, function(x) {if (all(x)) return ('?') else y <- x; y[lvls=='-']<-TRUE; if (nTokens > 4 && all(y)) return ('+') else return (output(lvls[x]))})
### }