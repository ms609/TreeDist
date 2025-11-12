# TreeDist: Distances between Phylogenetic Trees

'TreeDist' is an R package that implements a suite of metrics that
quantify the topological distance between pairs of unweighted
phylogenetic trees. It also includes a simple "Shiny" application to
allow the visualization of distance-based tree spaces, and functions to
calculate the information content of trees and splits.

## Details

"TreeDist" primarily employs metrics in the category of "generalized
Robinson–Foulds distances": they are based on comparing splits
(bipartitions) between trees, and thus reflect the relationship data
within trees, with no reference to branch lengths. Detailed
documentation and usage instructions are [available
online](https://ms609.github.io/TreeDist/) or in the vignettes.

### Generalized RF distances

The [Robinson–Foulds
distance](https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html)
simply tallies the number of non-trivial splits (sometimes inaccurately
termed clades, nodes or edges) that occur in both trees – any splits
that are not perfectly identical contributes one point to the distance
score of zero, however similar or different they are. By overlooking
potential similarities between almost-identical splits, this
conservative approach has undesirable properties.

["Generalized" RF
metrics](https://ms609.github.io/TreeDist/articles/Generalized-RF.html)
generate *matchings* that pair each split in one tree with a similar
split in the other. Each pair of splits is assigned a similarity score;
the sum of these scores in the optimal matching then quantifies the
similarity between two trees.

Different ways of calculating the the similarity between a pair of
splits lead to different tree distance metrics, implemented in the
functions below:

- [[`MutualClusteringInfo()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)](https://ms609.github.io/TreeDist/reference/TreeDistance.html),
  [[`SharedPhylogeneticInfo()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)](https://ms609.github.io/TreeDist/reference/TreeDistance.html)

  - Smith (2020) scores matchings based on the amount of information
    that one partition contains about the other. The Mutual Phylogenetic
    Information assigns zero similarity to split pairs that cannot both
    exist on a single tree; The Mutual Clustering Information metric is
    more forgiving, and exhibits more desirable behaviour; it is the
    recommended metric for tree comparison. (Its complement,
    [[`ClusteringInfoDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)](https://ms609.github.io/TreeDist/reference/TreeDistance.html),
    returns a tree distance.)

- [[`NyeSimilarity()`](https://ms609.github.io/TreeDist/dev/reference/NyeSimilarity.md)](https://ms609.github.io/TreeDist/reference/NyeSimilarity.html)

  - Nye *et al.* (2006) score matchings according to the size of the
    largest split that is consistent with both of them, normalized
    against the Jaccard index. This approach is extended by Böcker *et
    al*. (2013) with the Jaccard–Robinson–Foulds metric (function
    [[`JaccardRobinsonFoulds()`](https://ms609.github.io/TreeDist/dev/reference/JaccardRobinsonFoulds.md)](https://ms609.github.io/TreeDist/reference/JaccardRobinsonFoulds.html)).

- [[`MatchingSplitDistance()`](https://ms609.github.io/TreeDist/dev/reference/MatchingSplitDistance.md)](https://ms609.github.io/TreeDist/reference/MatchingSplitDistance.html)

  - Bogdanowicz and Giaro (2012) and Lin *et al.* (2012) independently
    proposed counting the number of "mismatched" leaves in a pair of
    splits.
    [[`MatchingSplitInfoDistance()`](https://ms609.github.io/TreeDist/dev/reference/TreeDistance.md)](https://ms609.github.io/TreeDist/reference/TreeDistance.html)
    provides an information-based equivalent (Smith 2020).

The package also implements the variation of the path distance proposed
by Kendal and Colijn (2016) (function
[[`KendallColijn()`](https://ms609.github.io/TreeDist/dev/reference/KendallColijn.md)](https://ms609.github.io/TreeDist/reference/KendallColijn.html)),
approximations of the Nearest-Neighbour Interchange (NNI) distance
(function
[[`NNIDist()`](https://ms609.github.io/TreeDist/dev/reference/NNIDist.md)](https://ms609.github.io/TreeDist/reference/NNIDist.html);
following Li *et al.* (1996)), and calculates the size (function
[[`MASTSize()`](https://ms609.github.io/TreeDist/dev/reference/MASTSize.md)](https://ms609.github.io/TreeDist/reference/MASTSize.html))
and information content (function
[[`MASTInfo()`](https://ms609.github.io/TreeDist/dev/reference/MASTSize.md)](https://ms609.github.io/TreeDist/reference/MASTSize.html))
of the Maximum Agreement Subtree.

For an implementation of the Tree Bisection and Reconnection (TBR)
distance, see the package
'[TBRDist](https://ms609.github.io/TBRDist/index.html)'.

## Tree space analysis

Map tree spaces and readily visualize mapped landscapes, avoiding common
analytical pitfalls (Smith, forthcoming), using the inbuilt graphical
user interface:

    TreeDist::MapTrees()

Serious analysts should consult the
[vignette](https://ms609.github.io/TreeDist/articles/treespace.html) for
a command-line interface.

## References

- Böcker S, Canzar S, Klau GW (2013). “The generalized Robinson-Foulds
  metric.” In Darling A, Stoye J (eds.), *Algorithms in Bioinformatics.
  WABI 2013. Lecture Notes in Computer Science, vol 8126*, 156–169.
  Springer, Berlin, Heidelberg.
  [doi:10.1007/978-3-642-40453-5_13](https://doi.org/10.1007/978-3-642-40453-5_13)
  .

- Bogdanowicz D, Giaro K (2012). “Matching split distance for unrooted
  binary phylogenetic trees.” *IEEE/ACM Transactions on Computational
  Biology and Bioinformatics*, **9**(1), 150–160.
  [doi:10.1109/TCBB.2011.48](https://doi.org/10.1109/TCBB.2011.48) .

- Kendall M, Colijn C (2016). “Mapping phylogenetic trees to reveal
  distinct patterns of evolution.” *Molecular Biology and Evolution*,
  **33**(10), 2735–2743.
  [doi:10.1093/molbev/msw124](https://doi.org/10.1093/molbev/msw124) .

- Li M, Tromp J, Zhang L (1996). “Some notes on the nearest neighbour
  interchange distance.” In Goos G, Hartmanis J, Leeuwen J, Cai J, Wong
  CK (eds.), *Computing and Combinatorics*, volume 1090, 343–351.
  Springer, Berlin, Heidelberg. ISBN 978-3-540-61332-9
  978-3-540-68461-9,
  [doi:10.1007/3-540-61332-3_168](https://doi.org/10.1007/3-540-61332-3_168)
  .

- Lin Y, Rajan V, Moret BME (2012). “A metric for phylogenetic trees
  based on matching.” *IEEE/ACM Transactions on Computational Biology
  and Bioinformatics*, **4**(9), 1014–1022.
  [doi:10.1109/TCBB.2011.157](https://doi.org/10.1109/TCBB.2011.157) .

- Nye TMW, Liò P, Gilks WR (2006). “A novel algorithm and web-based tool
  for comparing two alternative phylogenetic trees.” *Bioinformatics*,
  **22**(1), 117–119.
  [doi:10.1093/bioinformatics/bti720](https://doi.org/10.1093/bioinformatics/bti720)
  .

- Smith MR (2020). “Information theoretic Generalized Robinson-Foulds
  metrics for comparing phylogenetic trees.” *Bioinformatics*,
  **36**(20), 5007–5013.
  [doi:10.1093/bioinformatics/btaa614](https://doi.org/10.1093/bioinformatics/btaa614)
  .

- Smith MR (2022). “Robust analysis of phylogenetic tree space.”
  *Systematic Biology*, **71**(5), 1255–1270.
  [doi:10.1093/sysbio/syab100](https://doi.org/10.1093/sysbio/syab100) .

## See also

Further documentation is available in the [package
vignettes](https://ms609.github.io/TreeDist/articles/), visible from R
using `vignette(package = "TreeDist")`.

Other R packages implementing tree distance functions include:

- [ape](https://emmanuelparadis.github.io/):

  - [`cophenetic.phylo()`](https://rdrr.io/pkg/ape/man/cophenetic.phylo.html):
    Cophenetic distance

  - [`dist.topo()`](https://rdrr.io/pkg/ape/man/dist.topo.html): Path
    (topological) distance, Robinson–Foulds distance.

- [phangorn](https://cran.r-project.org/package=phangorn)

  - `treedist()`: Path, Robinson–Foulds and approximate SPR distances.

- [Quartet](https://ms609.github.io/Quartet/): Triplet and Quartet
  distances, using the tqDist algorithm.

- [TBRDist](https://ms609.github.io/TBRDist/): TBR and SPR distances on
  unrooted trees, using the 'uspr' C library.

- [distory](https://cran.r-project.org/package=distory) (unmaintained):
  Geodesic distance

## Author

**Maintainer**: Martin R. Smith <martin.smith@durham.ac.uk>
([ORCID](https://orcid.org/0000-0001-5660-1727)) \[copyright holder,
programmer\]

Other contributors:

- Roy Jonker <roy_jonker@magiclogic.com> (LAP algorithm) \[programmer,
  copyright holder\]

- Yong Yang <yongyanglink@gmail.com> (LAP algorithm) \[contributor,
  copyright holder\]

- Yi Cao (LAP algorithm) \[contributor, copyright holder\]
