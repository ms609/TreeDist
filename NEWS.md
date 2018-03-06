# TreeSearch 0.1.0

## New functions
- Added NJTree function as shortcut to generate Neighbour-Joining tree from a dataset
- Add functions to allow recovery of all trees one rearrangement from that input

## Efficiency gains
- Separate out NNISwap functions to allow more efficient rearrangement of edgeLists
- [9002] Improve efficiency by using three-pass algorithm in place of four-pass precursor.

## Bug fixes
- [9003] User now able to specify value of concavity constant (was overriden to k = 4)
- [9003] Bootstrap replicates now scored correctly (and without warning) under implied weights

# TreeSearch 0.0.7

## Inapplicables:
- Integrated with this package (previously in `inapplicable`)
- Handle inapplicable data via API to Martin Brazeau's Morphy Phylogenetic Library

## Profile Parsimony:
- Integrated with this package (previously in `ProfileParsimony`)
- Faster calculation of concavity profiles in C
- Persistent memoization with R.cache

# TreeSearch 0.0.6
- First CRAN submission
