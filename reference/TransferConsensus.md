# Consensus tree minimizing transfer distance

Construct a consensus tree that minimizes the sum of transfer distances
to a set of input trees, using a greedy add-and-prune heuristic.

## Usage

``` r
TransferConsensus(
  trees,
  scale = TRUE,
  greedy = c("best", "first"),
  init = c("empty", "majority")
)
```

## Arguments

- trees:

  An object of class `multiPhylo`: the input trees. All trees must share
  the same tip labels.

- scale:

  Logical; if `TRUE` (default), use the scaled transfer distance
  (normalized by light-side size minus one). If `FALSE`, use the
  unscaled (raw Hamming) transfer distance.

- greedy:

  Character string specifying the greedy strategy: `"best"` (default)
  picks the single highest-benefit action at each step; `"first"` picks
  the first improving action encountered (faster but potentially lower
  quality).

- init:

  Character string specifying the initial consensus: `"empty"` (default)
  starts with no splits (purely additive); `"majority"` starts with the
  majority-rule consensus and refines.

## Value

A tree of class `phylo`.

## Details

Unlike the majority-rule consensus, which minimizes Robinson-Foulds
distance and can be highly unresolved when phylogenetic signal is low,
`TransferConsensus()` uses the finer-grained transfer distance (Lemoine
et al. 2018) to construct a more resolved consensus tree.

The algorithm pools all splits observed across input trees, computes
pairwise transfer distances between them, and greedily adds or removes
splits to minimize total transfer dissimilarity cost. The approach
follows Takazawa et al. (2026) , reimplemented for 'TreeDist'
infrastructure.

## References

Lemoine F, Domelevo Entfellner J, Wilkinson E, Correia D, Dávila Felipe
M, De Oliveira T, Gascuel O (2018). “Renewing Felsenstein's phylogenetic
bootstrap in the era of big data.” *Nature*, **556**(7702), 452–456.
[doi:10.1038/s41586-018-0043-0](https://doi.org/10.1038/s41586-018-0043-0)
.  
  
Takazawa Y, Takeda A, Hayamizu M, Gascuel O (2026). “Outperforming the
majority-rule consensus tree using fine-grained dissimilarity measures.”
*bioRxiv*.
[doi:10.64898/2026.03.16.712085](https://doi.org/10.64898/2026.03.16.712085)
.

## See also

Other summary trees:
[`MCITree()`](https://ms609.github.io/TreeDist/reference/MCITree.md)

## Examples
