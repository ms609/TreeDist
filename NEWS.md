# TreeSearch 0.1.0
- Added NJTree function as shortcut to generate Neighbour-Joining tree from a dataset
- Separate out NNISwap functions to allow more efficient rearrangement of edgeLists
- Add functions to allow recovery of all trees one rearrangement from that input
- [9002] Improve efficiency by using three-pass algorithm in place of four-pass precursor.

# TreeSearch 0.0.7

## Inapplicables:
- Integrated with this package (previously in `inapplicable`)
- Handle inapplicable data via API to Martin Brazeau's Morphy Phylogenetic Library

## Profile Parsimony:
- Integrated with this package (previously in `ProfileParsimony`)
- Faster calculation of concavity profiles in C
- Persistent memoization with R.cache

# TreeSearch 0.0.6
- Submitted to CRAN
