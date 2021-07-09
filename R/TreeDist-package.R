#' TreeDist: Distances between Phylogenetic Trees
#' 
#' 'TreeDist' is an R package that implements a suite of metrics that quantify the
#' topological distance between pairs of unweighted phylogenetic trees.
#' It also includes a simple 'Shiny' application to allow the visualization of
#' distance-based tree spaces, and functions to calculate the information content
#' of trees and splits.
#' 
#' 'TreeDist' primarily employs metrics in the category of
#' 'generalized Robinson&ndash;Foulds distances': they are based on comparing splits
#' (bipartitions) between trees, and thus reflect the relationship data within 
#' trees, with no reference to branch lengths.
#' Detailed documentation and usage instructions are 
#' [available online](https://ms609.github.io/TreeDist/) or in the vignettes.
#' 
#' 
#' ## Generalized RF distances
#' 
#' The [Robinson&ndash;Foulds distance](https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html)
#' simply tallies the number of non-trivial splits (sometimes inaccurately
#' termed clades, nodes or edges) that occur in both trees -- any splits that are
#' not perfectly identical contributes one point to the distance score of zero, 
#' however similar or different they are.
#' By overlooking potential similarities between almost-identical splits, 
#' this conservative approach has undesirable properties.
#' 
#' ['Generalized' RF metrics](https://ms609.github.io/TreeDist/articles/Generalized-RF.html)
#' generate _matchings_ that pair each split in one tree with a similar split 
#' in the other.
#' Each pair of splits is assigned a similarity score; the sum of these scores in
#' the optimal matching then quantifies the similarity between two trees.
#' 
#' Different ways of calculating the the similarity between a pair of splits
#' lead to different tree distance metrics, implemented in the functions below:
#'   
#' * [`MutualClusteringInfo()`](https://ms609.github.io/TreeDist/reference/TreeDistance.html), [`SharedPhylogeneticInfo()`](https://ms609.github.io/TreeDist/reference/TreeDistance.html)
#' 
#'   + Smith (2020) scores matchings based on the amount of information
#'     that one partition contains about the other.  The Mutual Phylogenetic
#'     Information assigns zero similarity to split pairs that cannot
#'     both exist on a single tree; The Mutual 
#'     Clustering Information metric is more forgiving, and exhibits more 
#'     desirable behaviour; it is the recommended metric for tree comparison.
#'     (Its complement, [`ClusteringInfoDistance()`](https://ms609.github.io/TreeDist/reference/TreeDistance.html), returns a tree 
#'       distance.)
#' 
#' * [`NyeSimilarity()`](https://ms609.github.io/TreeDist/reference/NyeSimilarity.html)
#' 
#'   + Nye _et al._ (2006) score matchings according to the size of the largest 
#'     split that is consistent with both of them, normalized against 
#'     the Jaccard index.  This approach is extended by B&ouml;cker _et al_. (2013)
#'     with the Jaccard&ndash;Robinson&ndash;Foulds metric (function 
#'     [`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/reference/JaccardRobinsonFoulds.html)).
#' 
#' * [`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/reference/MatchingSplitDistance.html)
#' 
#'   + Bogdanowicz and Giaro (2012) and  Lin _et al._ (2012) independently proposed
#'     counting the number of 'mismatched' leaves in a pair of splits.
#'     [`MatchingSplitInfoDistance()`](https://ms609.github.io/TreeDist/reference/TreeDistance.html)
#'     provides an information-based equivalent (Smith 2020).
#' 
#' 
#' The package also implements the variation of the path distance 
#' proposed by Kendal and Colijn (2016) (function
#' [`KendallColijn()`](https://ms609.github.io/TreeDist/reference/KendallColijn.html)),
#' approximations of the Nearest-Neighbour Interchange (NNI) distance (function
#' [`NNIDist()`](https://ms609.github.io/TreeDist/reference/NNIDist.html); 
#' following Li _et al._ (1996)), and calculates the size (function
#' [`MASTSize()`](https://ms609.github.io/TreeDist/reference/MASTSize.html)) and 
#' information content (function
#' [`MASTInfo()`](https://ms609.github.io/TreeDist/reference/MASTSize.html)) of the 
#' Maximum Agreement Subtree.
#' 
#' For an implementation of the Tree Bisection and Reconnection (TBR) distance, see 
#' the package '[TBRDist](https://ms609.github.io/TBRDist/index.html)'.
#' 
#' 
#' # Tree space analysis
#' 
#' Map tree spaces and readily visualize mapped landscapes, avoiding
#' common analytical pitfalls (Smith, forthcoming),
#' using the inbuilt graphical user interface:
#'   
#'   ```r
#' TreeDist::MapTrees()
#' ```
#' 
#' Serious analysts should consult the
#' [vignette](https://ms609.github.io/TreeDist/articles/treespace.html)
#' for a command-line interface.
#' 
#' 
#' @seealso
#' 
#' Further documentation is available in the 
#' [package vignettes](https://ms609.github.io/TreeDist/articles/), visible from
#' R using `vignette(package = 'TreeDist')`.
#' 
#' Other R packages implementing tree distance functions include:
#'   
#'   * [ape](http://ape-package.ird.fr/):
#'     - `cophenetic.phylo()`: Cophenetic distance
#'     - `dist.topo()`: Path (topological) distance, Robinson&ndash;Foulds distance.
#'   * [phangorn](https://cran.r-project.org/package=phangorn)
#'     - `treedist()`: Path, Robinson&ndash;Foulds and approximate SPR distances.
#'   * [Quartet](https://ms609.github.io/Quartet/): Triplet and Quartet distances, 
#'   using the tqDist algorithm.
#'   * [TBRDist](https://ms609.github.io/TBRDist/): TBR and SPR distances on 
#'   unrooted trees, using the 'uspr; C library.
#'   * [distory](https://cran.r-project.org/package=distory) (unmaintained): 
#'     Geodesic distance
#' 
#' @references
#' 
#' - \insertRef{Bocker2013}{TreeDist}
#' 
#' - \insertRef{Bogdanowicz2012}{TreeDist}
#' 
#' - \insertRef{Kendall2016}{TreeDist}
#' 
#' - \insertRef{Li1996}{TreeDist}
#' 
#' - \insertRef{Lin2012}{TreeDist}
#' 
#' - \insertRef{Nye2006}{TreeDist}
#' 
#' - \insertRef{SmithDist}{TreeDist}
#' 
#' - \insertRef{SmithSpace}{TreeDist}
#' 
#' @encoding UTF-8
#' @keywords internal
"_PACKAGE"

# Suppress "NOTE: Nothing imported from Rdpack":
#' @importFrom Rdpack reprompt
NULL

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
