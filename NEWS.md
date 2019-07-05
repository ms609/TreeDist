# TreeSearch 0.3.2.9004 (develoment)

## New features
 - `ClusteringInfo` Adopts Meila's (2007) clustering information as a tree distance measure (but is not yet production-ready)
 - Implement another information theoretic tree distance measure (Smith, in prep)
 
## Bug fixes
 - Fix miscalculation of `MatchingSplitDifference`
 - Fix miscalculation of `JointInformation`

## Enhancements
 - `ReadTntTree` supports named taxa (with `taxname=`)
 - Improve speed of tests (by increasing probability of false positives)


# TreeSearch 0.3.2

## Enhancements
 - Improve text, content and build speed of vignettes


# TreeSearch 0.3.1

## New features
 - `NyeTreeSimilarity` function implements the tree similarity metric of
   Nye _et al._ (2006)
 - `MatchingSplitDistance` function implementing the Matching Split distance of 
   Bogdanowicz & Giaro (2012)

## Bug fixes
 - Check whether input tree is bifurcating before attempting rearrangements,
   to avoid crashes on unsupported input

# TreeSearch 0.3.0

## New features
 - Implement an information theoretic tree distance measure (Smith, in prep)
 - Prepare for new random number generator in R3.6.0

## Deprecations
 - Function `TreeSplits` is deprecated; use `Tree2Splits` instead

## Bug fixes
 - Correct some mistakes in the documentation


# TreeSearch 0.2.2

## Bug fixes 
 - Correct vignette titles


# TreeSearch 0.2.1

## New features
 - `CollapseNodes` and `CollapseEdges` allow the creation of polytomies
 - `Tree2Splits` lists the bipartition splits implied by a tree topology

## Enhancements
 - `SplitFrequency` now supports larger trees
 - Can specify tip labels directly to `ReadTntTree`, to avoid reliance on
   generative file

## Bug fixes
 - Export missing functions


# TreeSearch 0.2.0

## New features
 - `RootTree` function allows rooting of tree on incompletely specified
    or single-taxon outgroup
 - `AllTBR` returns all trees one TBR rearrangement away
 - `TBRMoves` reports all possible TBR rearrangements
 - `Jackknife` conducts Jackknife resampling
 - `SplitFrequency` reports frequency of clades in a forest
 - `SupportColour` allows visual marking of support values
 - `ApeTime` reports the creation date of an ape-exported tree
 - `SortTree` flips nodes into a consistent left-right order
 - `AsBinary` supports 0
 
## Enhancements
 - [IW]RatchetConsensus renamed to [IW]MultiRatchet, giving a better description 
     of the function's purpose
 - Don't warn about missing EOL when reading Nexus or TNT files
 - Add new 12-colour colourblind-friendly palette
 - FitchSteps now supports datasets with tips not found in tree
 - Improve portability of function `ReadTntTree`

## Bug fixes
 - [IW]MultiRatchet now considers trees identical even if they've been hit 
   a different number of times


# TreeSearch 0.1.2

## Bug fixes
- Update MorphyLib library to fix C warnings
- Remove non-ASCII characters from data
- Disable slow-building and problematic vignette
- Use local copy of citation style when building vignettes


# TreeSearch 0.1.0

## New features
- Helper functions to read Nexus and TNT data and trees
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
