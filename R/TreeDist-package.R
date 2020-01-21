#' TreeDist:  Distances between Phylogenetic Trees
#' 
#' `TreeDist` implements a suite of metrics that quantify the topological 
#' distance between pairs of unweighted phylogenetic trees.
#' The metrics generally fall in the category of "generalized Robinson-Foulds
#' distances": they are based on comparing partitions between trees, and thus
#' reflect the relationship data within trees, with no reference to branch 
#' lengths.
#' Detailed documentationand usage instructions are 
#' [available online](https://ms609.github.org/TreeDist) or in the vignettes.
#' 
#' The [Robinson-Foulds distance](https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html)
#' simply tallies the number of splits (loosely, clades) that occur in both trees -- 
#'   any splits that are not perfectly identical are assigned a score of zero, however similar 
#' or different they are.
#' By overlooking potential similarities between almost-identical splits, 
#' this conservative approach has undesirable properties.
#' 
#' ['Generalized' RF metrics](https://ms609.github.io/TreeDist/articles/Generalized-RF.html)
#' pair each split in one tree with a similar split in the other.
#' Each pair of splits is assigned a similarity 
#' score; the sum of these scores in the optimal partition matching 
#' then describes the similarity between two trees. 
#' 
#' Different ways of calculating the the similarity between a pair of partitions 
#' lead to different tree distance metrics, implemented in the functions below:
#'   
#'   
#' * [`MutualClusteringInfo`](https://ms609.github.io/TreeDist/reference/TreeDistance.html), [`SharedPhylogeneticInfo`](https://ms609.github.io/TreeDist/reference/TreeDistance.html)
#' 
#'   + Smith (forthcoming) scores matchings based on the amount of information
#'     that one partition contains about the other.  The Mutual Phylogenetic
#'     Information imposes arboreal matching: i.e. a pair of splits that cannot
#'     both exist on a single tree are assigned zero similarity.  The Mutual 
#'     Clustering Information metric is more forgiving, and exhibits more 
#'     desirable behaviour; it is the recommended metric for tree comparison.
#'     (The corresponding
#'     [`ClusteringInfoDistance`](https://ms609.github.io/TreeDist/reference/TreeDistance.html)
#'     returns a tree distance.)
#' 
#' * [`NyeTreeSimilarity`](https://ms609.github.io/TreeDist/reference/NyeTreeSimilarity.html)
#' 
#'   + Nye _et al._ (2006) score matchings according to the size of the largest 
#'     split that is consistent with both of them, normalized against 
#'     the Jaccard index.  This approach is extended by B&ouml;cker _et al_. (2013)
#'     with the Jaccard-Robinson-Foulds metric (function 
#'     [`JaccardRobinsonFoulds`](https://ms609.github.io/TreeDist/reference/JaccardRobinsonFoulds.html)).
#' 
#' * [`MatchingSplitDistance`](https://ms609.github.io/TreeDist/reference/MatchingSplitDistance.html)
#' 
#'  + Bogdanowicz and Giaro (2012), and independently Lin _et al._ (2012), count 
#'     the number of 'mismatched' terminals in a pair of splits. An 
#'     information-based equivalent (Smith, forthcoming) is provided in the function
#'     [`MatchingSplitInfoDistance`](https://ms609.github.io/TreeDist/reference/TreeDistance.html).
#'     
#' The package also implements the variation of the path distance 
#' proposed by Kendal and Colijn (2016) (function 
#' [`KendallColijn`](https://ms609.github.io/TreeDist/reference/KendallColijn.html)),
#' approximations of the Nearest-Neighbour Interchange (NNI) distance (function
#' [`NNIDist`](https://ms609.github.io/TreeDist/reference/NNIDist.html); following
#' Li _et al._ (1996)), and calculates the size (function
#' [`MASTSize`](https://ms609.github.io/TreeDist/reference/MASTSize.html)) and 
#' information content (function
#' [`MASTInfo`](https://ms609.github.io/TreeDist/reference/MASTSize.html)) of the 
#' Maximum Agreement Subtree.
#' 
#' For an implementation of the Tree Bisection and Reconnection (TBR) distance, 
#' see [`TBRDist`](https://ms609.github.io/TBRDist/index.html).
#' 
#' 
#' @seealso
#' 
#' Further documentation is available in the 
#' [package vignettes](https://ms609.github.io/TreeDist/articles/), visible from
#' R using `vignette(package='TreeDist')`.
#' 
#' Other tree distance functions are implemented in:
#'   
#'   * [ape](http://ape-package.ird.fr/):
#'     - `cophenetic.phylo`: Cophenetic distance
#'     - `dist.topo`: Path (topological) distance, Robinson-Foulds distance.
#' * [phangorn](https://cran.r-project.org/package=phangorn)
#'     - `treedist`: SPR, Robinson-Foulds and path distances.
#' * [Quartet](http://ms609.github.io/Quartet/): Triplet and Quartet distances, 
#' using the tqDist algorithm.
#' * [TBRDist](http://ms609.github.io/TBRDist/): TBR and SPR distances on 
#' unrooted trees, using the `uspr` C library.
#' * [distory](https://cran.r-project.org/package=distory) (unmaintained): Geodesic distance
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
#' - \insertRef{Nye2006}{TreeDist}
#' 
#' - \insertRef{SmithDist}{TreeDist}
#' 
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
