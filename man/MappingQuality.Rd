% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trustworthiness.R
\name{MappingQuality}
\alias{MappingQuality}
\alias{ProjectionQuality}
\title{Quality of mapped distances}
\usage{
MappingQuality(original, mapped, neighbours = 10L)

ProjectionQuality(original, mapped, neighbours = 10L)
}
\arguments{
\item{original, mapped}{Square matrix or \code{dist} object containing
original / mapped pairwise distances.}

\item{neighbours}{Number of nearest neighbours to use in calculation.}
}
\value{
\code{MappingQuality()} returns a named vector of length four,
containing the entries: \code{Trustworthiness}, \code{Continuity}, \code{TxC}
(the product of these values), and \code{sqrtTxC} (its square root).
}
\description{
\code{MappingQuality()} calculates the trustworthiness and continuity
of mapped distances \insertCite{Venna2001,Kaski2003}{TreeDist}.
Trustworthiness measures, on a scale from 0--1,
the degree to which points that are nearby in a mapping are truly close
neighbours; continuity, the extent to which points that are truly nearby
retain their close spatial proximity in a mapping.
}
\examples{
library('TreeTools', quietly = TRUE, warn.conflict = FALSE)
trees <- as.phylo(0:10, nTip = 8)
distances <- ClusteringInfoDistance(trees)
mapping <- cmdscale(distances)
MappingQuality(distances, dist(mapping), 4)
}
\references{
\insertAllCited{}
}
\seealso{
Other tree space functions: 
\code{\link{MapTrees}()},
\code{\link{SpectralEigens}()},
\code{\link{median.multiPhylo}()}
}
\author{
Wrapper for functions from Charlotte Soneson's \pkg{dreval},
https://github.com/csoneson/dreval/blob/master/R/trustworthiness.R
}
\concept{tree space functions}
