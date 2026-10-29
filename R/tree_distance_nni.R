#' Approximate Nearest Neighbour Interchange distance
#' 
#' Use the approach of \insertCite{Li1996;textual}{TreeDist} to approximate the
#' Nearest Neighbour Interchange distance \insertCite{Robinson1971}{TreeDist} 
#' between phylogenetic trees.
#' 
#' In brief, this approximation algorithm works by identifying edges in one
#' tree that do not match edges in the second.  Each of these edges must
#' undergo at least one NNI operation in order to reconcile the trees.
#' Edges that match in both trees need never undergo an NNI operation, and 
#' divide each tree into smaller regions.  By "cutting" matched edges into two,
#' a tree can be divided into a number of regions that solely comprise unmatched
#' edges.
#' 
#' These regions can be viewed as separate trees that need to be reconciled.
#' One way to reconcile these trees is to conduct a series of NNI operations
#' that reduce a tree to a pectinate (caterpillar) tree, then to conduct an 
#' analogue of the mergesort algorithm.  This takes at most _n_ log _n_ + O(_n_)
#' NNI operations, and provides a loose upper bound on the NNI score. 
#' The maximum number of moves for an _n_-leaf tree
#' ([OEIS A182136](https://oeis.org/A182136)) can be calculated exactly for
#' small trees \insertCite{Fack2002}{TreeDist}; this provides a tighter upper
#' bound, but is unavailable for _n_ > 12.
#' `NNIDiameter()` reports the limits on this bound.
#' 
#' Where `regionTips` is positive, each region with up to `regionTips` leaves
#' instead contributes the exact NNI distance between its two subtrees,
#' found by bidirectional breadth-first search of tree space.
#' Reconciling each region separately is one way to transform one tree into
#' the other, so the sum, `region_upper`, is an upper bound on the distance.
#' It need not be the exact distance: a shortest path may break a split
#' that the two trees share \insertCite{Li1996}{TreeDist}, so NNI distance
#' does not decompose into regions.
#' Larger regions, and any region whose search would visit more than two
#' million trees, contribute their `best_upper` value.
#' Searching a region of up to 11 leaves takes under a second; larger regions
#' repay the search only when their two subtrees are similar, as a search
#' that exhausts its budget takes several seconds.
#' 
#' 
#' \tabular{rccccccccccccc}{
#'   Leaves:   \tab 1 \tab 2 \tab 3 \tab 4 \tab 5 \tab 6 \tab 7 \tab 8 \tab 9
#'    \tab 10 \tab 11 \tab 12 \tab 13 \cr
#'   Diameter: \tab 0 \tab 0 \tab 0 \tab 1 \tab 3 \tab 5 \tab 7 \tab 10 \tab 12
#'    \tab 15 \tab 18 \tab 21 \tab ? \cr
#'  }
#' 
#' @template tree12Params
#' @param regionTips Integer specifying the largest region, in leaves, whose
#' distance is computed exactly when calculating `region_upper`, up to 16;
#' `0` omits `region_upper`.
#' 
#' @return `NNIDist()` returns, for each pair of trees, a named vector
#'  containing six integers, or seven if `regionTips > 0`:
#'
#' - `lower` is a lower bound on the NNI distance: the number of edges in
#' one tree that do not occur in the other, i.e. half the Robinson-Foulds
#' distance between the trees.
#'
#' - `best_upper` is the tightest available upper bound: it uses the exact
#' diameter for regions with < 13 leaves, and `loose_upper` for larger
#' regions.
#'
#' - `tight_upper` is an upper bound on the distance, based on calculated
#' maximum diameters for trees with < 13 leaves.  _NA_ is returned if trees are
#' too different to employ this approach.
#'
#' - `loose_upper` is a looser upper bound on the distance, using _n_ log _n_ +
#' O(_n_): the smaller of `fack_upper` and `li_upper`.
#'
#' - `fack_upper` is an upper bound based on the maximum diameter given by
#' theorem 3 of \insertCite{Fack2002;textual}{TreeDist}.
#'
#' - `li_upper` is an upper bound based on the sorting approach of
#' \insertCite{Li1996;textual}{TreeDist}.
#'
#' - `region_upper` is an upper bound that sums the exact distances between
#' the two trees' subtrees within each region of up to `regionTips` leaves
#' (see Details); it is never larger than `best_upper`.
#'
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @examples
#' library("TreeTools", quietly = TRUE)
#' 
#' NNIDist(BalancedTree(7), PectinateTree(7))
#' 
#' NNIDist(BalancedTree(7), as.phylo(0:2, 7))
#' NNIDist(as.phylo(0:2, 7), PectinateTree(7))
#'
#' NNIDist(list(bal = BalancedTree(7), pec = PectinateTree(7)),
#'         as.phylo(0:2, 7))
#'
#' CompareAll(as.phylo(30:33, 8), NNIDist)
#'
#' NNIDist(BalancedTree(9), PectinateTree(9), regionTips = 9)
#' @template MRS
#'   
#' @family tree distances
#' @export
NNIDist <- function(tree1, tree2 = tree1, regionTips = 0L) {
  .TreeDistance(.NNIDistSingle, tree1, tree2, regionTips = regionTips)
}

#' @importFrom TreeTools Postorder RenumberTips
#' @importFrom ape Nnode.phylo
.NNIDistSingle <- function(tree1, tree2, nTip, regionTips = 0L, ...) {
  .CheckMaxTips(nTip)
  if (nrow(tree1[["edge"]]) != nrow(tree2[["edge"]])) {
    stop("Both trees must have the same number of edges. ",
         "Is one rooted and the other unrooted?")
  }
  tree2 <- RenumberTips(tree2, as.character(tree1[["tip.label"]]))
  
  edge1 <- Postorder(tree1[["edge"]])
  edge2 <- Postorder(tree2[["edge"]])
  
  cpp_nni_distance(edge1, edge2, nTip, regionTips)
}

# Exact NNI distance between two binary trees of up to 16 leaves, by
# breadth-first search; NA if more than `budget` trees must be visited.
#' @importFrom TreeTools Postorder RenumberTips NTip
.NNIExact <- function(tree1, tree2, budget = 2^22) {
  tree2 <- RenumberTips(tree2, tree1[["tip.label"]])
  cpp_nni_exact(
    Postorder(tree1)[["edge"]],
    Postorder(tree2)[["edge"]],
    NTip(tree1),
    budget
  )
}

#' Bounds on the diameter of the NNI distance
#' 
#' @param tree Object of supported class representing a tree or list of trees,
#' or an integer specifying the number of leaves in a tree/trees.
#' 
#' @return `NNIDiameter()` returns a matrix specifying (bounds on) the diameter
#' of the NNI distance metric on the specified tree(s).
#' Columns correspond to:
#' 
#' - `liMin`:  \deqn{n - 3}{_n_ &minus; 3;}, a lower bound on the diameter
#'   (Li _et al._ 1996);
#'   
#' - `fackMin`: Lower bound on diameter following Fack _et al_. (2002), i.e.
#'   \deqn{\log2{N!} / 4}{log&#8322; (_N_!) / 4};
#'   
#' - `min`: The larger of `liMin` and `fackMin`;
#' 
#' - `exact`: The exact value of the diameter, where _n_ &lt; 13;
#' 
#' - `liMax`: Upper bound on diameter following Li _et al._ (1996), i.e. 
#'   \deqn{n \log2{n} + \textrm{O}(n)}{n log&#8322; _n_ + O(_n_)};
#'   
#' - `fackMax`: Upper bound on diameter following theorem 3 of
#'   Fack _et al_. (2002), i.e. _N_ _c_ &minus; 2^_c_ + 1 + 2(_N_ &minus;
#'   ceiling(log&#8322; (_N_ + 1))), where _c_ = ceiling(log&#8322; _N_);
#'   
#' - `max`: The smaller of `liMax` and `fackMax`;
#'   
#' where _n_ is the number of leaves, and _N_ the number of internal nodes,
#' i.e. \deqn{n - 2}{_n_ &minus; 2}.
#'   
#' @encoding UTF-8
#' @rdname NNIDist
#' @export
NNIDiameter <- function(tree) UseMethod("NNIDiameter")

.SortingNumber <- function(n_tip) {
  lgN <- ceiling(log2(n_tip))
  n_tip * lgN - 2^lgN + 1
}

.DegenerateDistance <- function(n_tip) { # For working, see nni_distance.cpp
  nodes_in_full <- pmax(ceiling(log2(n_tip) - log2(3)), 0L)
  tips_left <- n_tip - 2^nodes_in_full
  min_backbone_nodes <- pmax(nodes_in_full + ceiling(log2(tips_left)), 0L)
  
  # Return:
  pmax(n_tip - 2, 0) - min_backbone_nodes
}

#' @export
NNIDiameter.numeric <- function(tree) {
  n <- tree - 2L # unrooted tree with n + 2 leaves
  n[tree < 3L] <- NA
  
  liMin <- tree - 3L
  fackMin <- ceiling(lfactorial(n) / 4L / log(2))
  
  liMax <- .SortingNumber(tree) + (2 * .DegenerateDistance(tree))
  lgN <- ceiling(log2(n))
  fackMax <- (n * lgN) - (2 ^ lgN) + 1 + (2 * (n - ceiling(log2(n + 1))))
  cbind(liMin = liMin,
        fackMin = fackMin,
        min = pmax(liMin, fackMin),
        exact = c(0, 1, 3, 5, 7, 10, 12, 15, 18, 21)[n],
        liMax = liMax,
        fackMax = fackMax,
        max = pmin(liMax, fackMax)
  )
}

#' @importFrom TreeTools NTip
#' @export
NNIDiameter.phylo <- function(tree) {
  NNIDiameter(NTip(tree))
}

#' @export
NNIDiameter.multiPhylo <- NNIDiameter.phylo

#' @export
NNIDiameter.list <- function(tree) {
  lapply(tree, NNIDiameter)
}
