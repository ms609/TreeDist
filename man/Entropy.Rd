% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_distance_utilities.R
\name{Entropy}
\alias{Entropy}
\title{Entropy in bits}
\usage{
Entropy(...)
}
\arguments{
\item{\dots}{Numerics or numeric vector specifying probabilities of outcomes.}
}
\value{
\code{Entropy()} returns the entropy of the specified probabilities,
in bits.
}
\description{
Calculate the entropy of a vector of probabilities, in bits.
Probabilities should sum to one.
Probabilities equalling zero will be ignored.
}
\examples{
Entropy(1/2, 0, 1/2) # = 1
Entropy(rep(1/4, 4)) # = 2
}
\author{
\href{https://orcid.org/0000-0001-5660-1727}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
