# TreeSearch 0.1.1

## Bug fixes
- Update MorphyLib library to fix C warnings
- Remove non-ASCII characters from data


# TreeSearch 0.1.0

## New features
- Helper functions to read Nexus and TNT data and trees.
- Brewer palette in local data to allow easier colouring

## Enhancements
- Allow additional parameters to be passed to `consensus` via `ConsensusWithout`

## Bug fixes
- C11 compliance
- `IWRatchetConsensus` now relays concavity value to subsequent functions
- `ReadCharacters` returns labels for all characters and states if `character_num = NULL`


# TreeSearch 0.0.8

## New features
- Added NJTree function as shortcut to generate Neighbour-Joining tree from a dataset
- Add functions to allow recovery of all trees one rearrangement from that input

## Efficiency gains
- Separate out NNISwap functions to allow more efficient rearrangement of edgeLists
- [9002] Improve efficiency by using three-pass algorithm in place of four-pass precursor
- [9004] Bootstrap search improvements

## Bug fixes
- [9003] User now able to specify value of concavity constant (was overridden to k = 4)
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
