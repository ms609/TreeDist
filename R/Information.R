#' Mutual information of two splits
#' 
#' Reports the mutual phylogenetic information, or variation of phylogenetic
#' information, of two splits.
#' 
#' The mutual phylogenetic information corresponds to the entropy of the subset
#' of trees consistent with both splits; two splits that are consistent with a 
#' smaller number of trees will have a higher mutual information content.
#' 
#' Split 1 divides `n` terminals into two partitions, _A1_ and _B1_.
#' Split 2 divides the same terminals into the partitions _A2_ and _B2_.
#' 
#' Partitions must be named such that _A1_ fully overlaps with _A2_: 
#' that is to say, all taxa in _A1_ are also in _A2_, or _vice versa_.
#' Thus, all taxa in the smaller of _A1_ and _A2_ also occur in the larger.
#' 
#' @param n Integer specifying the number of terminals.
#' @param A1,A2 Integers specifying the number of taxa in _A1_ and _A2_, 
#' once the splits have been arranged such that _A1_ fully overlaps with _A2_.
#' @return 
#' `TreesConsistentWithTwoSplits` returns the number of unrooted bifurcating
#' trees consistent with two splits.
#' 
#' `SplitMutualInformation` returns the information that two splits 
#' have in common, in bits.
#' 
#' `SplitVariationOfInformation` returns the variation of information (Meila 2007)
#'  between the two splits, a measure of their difference, in bits.
#' 
#' @examples 
#'   # Eight tips, labelled A to H.
#'   # Split 1: ABCD:EFGH
#'   # Split 2: ABC:DEFGH
#'   # Let A1 = ABCD (four taxa), and A2 = ABC (three taxa).
#'   # A1 and A2 overlap (both contain ABC).
#'   
#'   TreesConsistentWithTwoSplits(n=8, A1=4, A2=3)
#'   SplitMutualInformation(n=8, A1=4, A2=3)
#'   SplitVariationOfInformation(n=8, A1=4, A2=3)
#'
#'   # If splits are identical, then their mutual information is the same
#'   # as the information of either split:
#'   SplitMutualInformation(n=8, A1=3, A2=3)
#'   TreeSearch::SplitInformation(3, 5)
#'   
#' @references \insertRef{Meila2007}{TreeDist}
#' 
#' @author Martin R. Smith
#' @family split information functions
#' @importFrom TreeSearch LogTreesMatchingSplit LnUnrooted
#' @export
SplitMutualInformation <- function(n, A1, A2 = A1) {
  (LogTreesMatchingSplit(A1, n - A1) 
   + LogTreesMatchingSplit(A2, n - A2)
   - LogTreesConsistentWithTwoSplits(n, A1, A2)
   - LnUnrooted(n)) / -log(2)
}

#' @describeIn SplitMutualInformation Variation of information between two splits.
#' @importFrom TreeSearch SplitInformation
#' @export
SplitVariationOfInformation <- function (n, A1, A2 = A1) {
  # TODO calculate more efficiently from first principles
  mutual <- SplitMutualInformation(n, A1, A2)
  SplitInformation (A1, n - A1) - mutual +
    SplitInformation(A2, n - A2) - mutual
}

#' Use Variation of Information to compare partitions as clusterings
#' 
#' Compare a pair of splits solely as clusterings of taxa, disregarding their
#' phlogenetic information, using the variation of clustering information 
#' proposed by Meila (2007).
#' 
#' This is equivalent to the mutual clustering information (Vinh 2010).
#' 
#' Whereas Meila and Vinh _et al_. are concerned with the entropy, i.e. the 
#' number of bits required to encode each tip, we are interested in the 
#' total information content, i.e. the entropy multiplied by the number of tips.
#' This is the value that the function returns.
#' 
#' @return Variation of information, measured in bits.
#' 
#' @references {
#'   \insertRef{Meila2007}{TreeDist}
#'   \insertRef{Vinh2010}{TreeDist}
#' }
#' 
#' @examples 
#'   # Maximum variation = information content of each split separately
#'   A <- TRUE
#'   B <- FALSE
#'   MeilaVariationOfInformation(c(A,A,A,B,B,B), c(A,A,A,A,A,A))
#'   Entropy(c(3, 3) / 6) + Entropy(c(0, 6) / 6)
#'   
#'   # Minimum variation = 0
#'   MeilaVariationOfInformation(c(A,A,A,B,B,B), c(A,A,A,B,B,B))
#'   
#'   # Not always possible for two evenly-sized splits to reach maximum
#'   # variation of information
#'   Entropy(c(3, 3) / 6) * 2  # = 2
#'   MeilaVariationOfInformation(c(A,A,A,B,B,B), c(A,B,A,B,A,B)) # < 2
#'   
#'   # Phylogenetically uninformative groupings contain partitioning information
#'   Entropy(c(1, 5) / 6)
#'   MeilaVariationOfInformation(c(B,A,A,A,A,A), c(A,A,A,A,A,B))
#' 
#' 
#' @template split12Params
#' @author Martin R. Smith
#' @export
MeilaVariationOfInformation <- function (split1, split2) {
  if (length(split1) != length(split2)) stop("Split lengths differ")
  n <- length(split1)
  
  associationMatrix <- c(sum(split1 & split2), sum(split1 & !split2),
                         sum(!split1 & split2), sum(!split1 & !split2))
  probabilities <- associationMatrix / n
  p1 <- probabilities[1] + probabilities[3]
  p2 <- probabilities[1] + probabilities[2]
  
  jointEntropies <- Entropy(probabilities)
  
  variationOfInformation <- jointEntropies + jointEntropies - 
    Entropy(c(p1, 1 - p1)) -
    Entropy(c(p2, 1 - p2))
  # Return:
  n * variationOfInformation
}
#' @describeIn MeilaVariationOfInformation Mutual clustering information of two splits
#' @return Mutual information, measured in bits.
#' @export
MeilaMutualInformation <- function (split1, split2) {
  if (length(split1) != length(split2)) stop("Split lengths differ")
  n <- length(split1)
  
  associationMatrix <- c(sum(split1 & split2), sum(split1 & !split2),
                         sum(!split1 & split2), sum(!split1 & !split2))
  probabilities <- associationMatrix / n
  p1 <- probabilities[1] + probabilities[3]
  p2 <- probabilities[1] + probabilities[2]
  
  jointEntropies <- Entropy(probabilities)
  mutualInformation <- Entropy(c(p1, 1 - p1)) + Entropy(c(p2, 1 - p2)) -
    jointEntropies
  
  # Return:
  if (abs(mutualInformation) < 1e-15) 0 else n * mutualInformation
}

# # #' Probability of matching this well
# # #' 
# # #' Probability that two random splits of the sizes provided will be at least
# # #' as similar as the two specified.
# # #'
# # #' @template split12Params
# # #' 
# # #' @return The proportion of permissable informative splits 
# # #' splitting the terminals into bipartitions of the sizes given,
# # #'  that match as well as `split1` and `split2` do.
# # #'  
# # AnySizeSplitMatchProbability <- function (split1, split2) {
# #   
# #   observedEntropy <- SplitEntropy(split1, split2)
# #   
# #   possibleEntropies <- AllEntropies(length(split1))
# #   
# #   nArrangements <- iMax - minA1B2 + 1L
# # 
# #   arrangements <- vapply(minA1B2:iMax,
# #                          function (i) c(A1 - i, i, i + A2 - A1, B2 - i),
# #                          double(4))
# #   
# #   H <- function(p) -sum(p[p > 0] * log(p[p > 0]))
# #   jointEntropies <- apply(arrangements / n, 2, H)
# #   # mutualInformation <- H(c(A1,n - A1) / n) + H(c(A2, B2) / n) - jointEntropies
# #   
# #   # Meila 2007; less is best
# #   variationOfInformation <- jointEntropies+ jointEntropies - 
# #     (H(c(A1,n - A1) / n) + H(c(A2, B2) / n))
# #   
# #   choices <- apply(arrangements, 2, NPartitionPairs)
# #   
# #   # Return:
# #   sum(choices[ranking <= ranking[partitions[1, 2] + 1L - minA1B2]]) / choose(n, A1)
# # }


#' Probability of matching this well
#' 
#' Calculates the probability that two random splits of the sizes provided
#' will be at least as similar as the two specified.
#'
#' @template split12Params
#' 
#' @return The proportion of permissable informative splits 
#' splitting the terminals into bipartitions of the sizes given,
#'  that match as well as `split1` and `split2` do.
#'  
#' @examples 
#' SplitMatchProbability(split1 = c(rep(TRUE, 4), rep(FALSE, 4)),
#'                       split2 = c(rep(TRUE, 3), rep(FALSE, 5)))
#'  
#' @author Martin R. Smith
#' @export
SplitMatchProbability <- function (split1, split2) {
  
  if (length(split1) != length(split2)) stop("Split lengths differ")
  
  partitions <- matrix(c(sum(split1 & split2),
                         sum(split1 & !split2),
                         sum(!split1 & split2),
                         sum(!split1 & !split2)
                         ), 2, 2) 
  #, dimnames=list(c('A1', 'B1'), c('A2', 'B2')))
  
  split1Size <- rowSums(partitions)
  split2Size <- colSums(partitions)
  A1 <- split1Size[1]
  A2 <- split2Size[1]
  B2 <- split2Size[2]
  n <- sum(partitions)
  
  minA1B2 <- max(0, A1 - A2)
  iMax <- min(A1, B2)
  minB1A2 <- B2 - iMax
  
  nArrangements <- iMax - minA1B2 + 1L
  # It turns out that this is called a confusion matrix, association matrix
  # or contingency table; see Meila 2007
  arrangements <- vapply(minA1B2:iMax,
                         function (i) c(A1 - i, i, i + A2 - A1, B2 - i),
                         double(4))
  
  #H <- function(p) -sum(p[p > 0] * log(p[p > 0]))
  #jointEntropies <- apply(arrangements / n, 2, H)
  #mutualInformation <- H(c(A1,n - A1) / n) + H(c(A2, B2) / n) - jointEntropies
  
  extraTipsInPerfectMatch.A1A2.B1B2 <- sum(arrangements[c(1, 4), nArrangements])
  extraTipsInPerfectMatch.A1B2.B1A2 <- sum(arrangements[c(2, 3), 1L])
  
  ranking <- if (extraTipsInPerfectMatch.A1A2.B1B2 > 
                 extraTipsInPerfectMatch.A1B2.B1A2) {
    c(seq.int(1L, nArrangements, 2L),
      rev(seq.int(2L, nArrangements, 2L)))
  } else if (extraTipsInPerfectMatch.A1A2.B1B2 < 
             extraTipsInPerfectMatch.A1B2.B1A2) {
    c(seq.int(2L, nArrangements, 2L),
      rev(seq.int(1L, nArrangements, 2L)))
  } else {
    c(seq.int(1L, nArrangements, 2L),
      rev(seq.int(2L, nArrangements, 2L)) - 1)
  }
  
  choices <- apply(arrangements, 2, NPartitionPairs)
  
  # Return:
  sum(choices[ranking <= ranking[partitions[1, 2] + 1L - minA1B2]]) /
    choose(n, A1)
}


#' All split pairings
#' 
#' Calculates the variation of (clustering) information (Meila 2007)
#' for each possible pairing of non-trivial splits,
#' and tabulates the number of pairings with each similarity.
#' 
#' Trivial splits -- those that divide one or zero taxa from the rest --
#' are not considered.
#' 
#' @param n Integer specifying the number of terminal taxa.
#' 
#' @return `AllSplitPairings` returns a named vector, specifying the number of 
#' split pairings producing
#' the variation of information given (in bits) in the name.  Splits
#' AB:CD and CD:AB are treated as distinct, so division of all values by four 
#' is justified in cases where unique pairings only are required.
#' 
#' @examples
#' AllSplitPairings(5)
#' 
#' @references 
#' * \insertRef{Meila2007}{TreeDist}
#' * \insertRef{SmithDist}{TreeDist}
#' 
#' @author Martin R. Smith
#' @importFrom memoise memoise
#' @export
AllSplitPairings <- memoise(function (n) {
  
  if (n < 4L) stop("No informative splits with < 4 taxa")
  
  # smallHalves <- 1L + seq_len(ceiling(n / 2) - 2L)
  dataRows <- 2L

  unevenPairs <- matrix(
    # For i in 2:largestSmallSplit
    # TODO: Make faster by not calculating bottom triangle
    unlist(lapply(1L + seq_len(n - 3L), function (inA) {
      # For j in 2:(n - 2)
      nCa <- choose(n, inA)
      outA <- n - inA
      hA <- Entropy(c(inA, outA) / n)
      unlist(lapply(1L + seq_len(n - 3L), function (inB) {
        outB <- n - inB
        hB <- Entropy(c(inB, outB) / n)
        vapply(max(0, inA + inB - n):min(inA, inB), function (inAB) {
          association <- c(inAB, inA - inAB, inB - inAB, n + inAB - inA - inB)
          jointEntropies <- Entropy(association / n)
          
          c(#inA, inB, inAB, 
            #npairs = NPartitionPairs(association), nis = choose(n, i),
            nTotal = nCa * choose(inA, inAB) * choose(outA, inB - inAB),
            VoI = jointEntropies + jointEntropies - hA - hB)
        }, double(dataRows))
      }))
    })), dataRows, dimnames=list(c('nTotal', 'VoI'), NULL))
  
  tapply(unevenPairs['nTotal', ], unevenPairs['VoI', ], sum)
})

#' Distributions of taxa consistent with a partition pair.
#' 
#' Number of terminal arrangements matching a specified configuration of 
#' two partitions.
#' 
#' Consider partitions that divide eight terminals, labelled A to H.
#' 
#' \tabular{rcll}{
#'   Bipartition 1:\tab ABCD:EFGH\tab A1 = ABCD\tab B1 = EFGH \cr
#'   Bipartition 2:\tab ABE:CDFGH\tab A2 = ABE\tab B2 = CDFGH
#' }
#' 
#' This can be represented by an association matrix:
#' 
#' \tabular{rll}{
#'      \tab *A2* \tab *B2* \cr
#' *A1* \tab AB  \tab C   \cr
#' *B1* \tab E   \tab FGH
#' }
#' 
#' The cells in this matrix contain 2, 1, 1 and 3 terminals respectively; this
#' four-element vector (`c(2, 1, 1, 3)`) is the `configuration` implied by 
#' this pair of bipartition splits.
#' 
#' @param configuration Integer vector of length four specifying the number of 
#' terminals that occur in both 
#'   (1) splits A1 and A2; 
#'   (2) splits A1 and B2;
#'   (3) splits B1 and A2;
#'   (4) splits B1 and B2.
#'   
#' @return The number of ways to distribute `sum(configuration)` taxa according
#'  to the specified pattern.
#'
#' @examples 
#' NPartitionPairs(c(2, 1, 1, 3))
#'  
#' @author Martin R. Smith
#' @export
NPartitionPairs <- function (configuration) {
  choose(sum(configuration[c(1, 3)]), configuration[1]) *
  choose(sum(configuration[c(2, 4)]), configuration[2])
}
#' @describeIn SplitMatchProbability The natural logarithm of the probability
LnSplitMatchProbability <- function(split1, split2) {
  log(SplitMatchProbability(split1, split2))
}

#' Entropy of two splits
#' 
#' Reports various values pertaining to the phylogenetic information content 
#' of two splits,
#' treating splits as subdivisions of _n_ terminals into two clusters.
#' 
#' @template split12Params
#' 
#' @references 
#' \insertRef{Meila2007}{TreeDist}
#' 
#' @return A numeric vector listing, in bits,
#'  * `h1` The entropy of split 1
#'  * `h1` The entropy of split 2
#'  * `jointH` The joint entropy of both splits
#'  * `i` The mutual information of the splits
#'  * `vI` The variation of information of the splits (see Meila 2007)
#' 
#' @author Martin R. Smith
#' @family split information functions
#' @export
SplitEntropy <- function (split1, split2=split1) {
  A1A2 <- sum(split1 & split2)
  A1B2 <- sum(split1 & !split2)
  B1A2 <- sum(!split1 & split2)
  B1B2 <- sum(!split1 & !split2)
  overlaps <- c(A1A2, A1B2, B1A2, B1B2)
  
  A1 <- A1A2 + A1B2
  A2 <- A1A2 + B1A2
  B1 <- B1A2 + B1B2
  B2 <- A1B2 + B1B2
  n <- A1 + B1
  
  h1 <- Entropy(c(A1, B1) / n)
  h2 <- Entropy(c(A2, B2) / n)
  jointH <- Entropy(overlaps[overlaps > 0L] / n)
  mutualInformation <- h1 + h2 - jointH
  variationOfInformation <- jointH - mutualInformation
  
  # Return:
  c(h1 = h1, h2 = h2, jointH = jointH, i = mutualInformation,
    vI = variationOfInformation)
}

#' Joint Information of two splits
#' 
#' Calculates the joint phylogenetic information content of two splits: 
#' that is, the information content of the two splits considered separately,
#' minus their mutual information (which would otherwise be counted twice).
#'
#' Because some information is common to both splits (`SplitMutualInformation`),
#' the joint information of two splits will be less than the sum of the
#' information of the splits taken separately -- unless the splits are
#' contradictory.
#' 
#' Split Y1 is defined as dividing taxa into the two sets A1 and B1,
#' and Y2=A2:B2.
#' 
#' Consider partitions that divide eight terminals, labelled A to H.
#' 
#' \tabular{rcll}{
#'   Bipartition 1:\tab ABCD:EFGH\tab A1 = ABCD\tab B1 = EFGH \cr
#'   Bipartition 2:\tab ABC:DEFGH\tab A2 = ABC\tab B2 = DEFGH
#' }
#' 
#' This can be represented by an association matrix:
#' 
#' \tabular{rll}{
#'       \tab *A2* \tab *B2* \cr
#'  *A1* \tab ABC  \tab D   \cr
#'  *B1* \tab      \tab EFGH
#' }
#' 
#' The joint information is given by
#' `JointInformation(3, 1, 0, 4)`.
#'
#' @param A1A2,A1B2,B1A2,B1B2 Number of taxa common to splits A1 and A2 (etc.).
#' 
#' @return `JointInformation` returns the joint phylogenetic information content
#' of two splits, measured in bits.
#' 
#' @examples 
#' JointInformation(3, 1, 0, 4)
#' 
#' @author Martin R. Smith
#' @family split information functions
#' @importFrom TreeSearch NUnrooted
#' @export
JointInformation <- function(A1A2, A1B2, B1A2, B1B2) {
  # Y1 = A1:B1
  A1 <- A1A2 + A1B2
  B1 <- B1A2 + B1B2
  n  <- A1 + B1
  
  # Y2 = A2:B2
  A2 <- A1A2 + B1A2
  B2 <- A1B2 + B1B2
  
  mutualInformation <- if(A1B2 == 0 || B1A2 == 0) {
    -log2(TreesConsistentWithTwoSplits(n, A1, A2) / NUnrooted(n))
  } else if(A1A2 == 0 || B1B2 == 0) {
    -log2(TreesConsistentWithTwoSplits(n, A1, B2) / NUnrooted(n))
  } else {
    0
  }
  
  # Return:
  SplitInformation(A1, B1) + SplitInformation(A2, B2) - mutualInformation
}

#' @describeIn SplitMutualInformation Number of trees consistent with two splits.
#' @family split information functions
#' @importFrom TreeSearch TreesMatchingSplit NRooted
#' @export
TreesConsistentWithTwoSplits <- function (n, A1, A2=A1) {
  
  smallSplit <- min(A1, A2)
  bigSplit <- max(A1, A2)
  
  if (smallSplit == 0) return (TreesMatchingSplit(bigSplit, n - bigSplit))
  if (bigSplit == n) return (TreesMatchingSplit(smallSplit, n - smallSplit))
  
  overlap <- bigSplit - smallSplit
  
  #  Here are two spits:
  #  AA OOO BBBBBB
  #  11 111 000000
  #  11 000 000000
  #  
  #  There are (2O - 5)!! unrooted trees of the overlapping taxa (O)
  #  There are (2A - 5)!! unrooted trees of A
  #  There are (2B - 5)!! unrooted trees of B
  #  
  #  There are 2A - 3 places on A that A can be attached to O.
  #  There are 2O - 3 places on O to which A can be attached.
  #  
  #  There are 2B - 3 places in B that B can be attached to (O+A).
  #  There are 2O - 3 + 2 places on O + A to which B can be attached:
  #   2O - 3 places on O, plus the two new edges created when A was joined to O.
  #  
  #  (2A - 3)(2A - 5)!! == (2A - 3)!!
  #  (2B - 3)(2B - 5)!! == (2B - 3)!!
  #  (2O - 1)(2O - 3)(2O - 5)!! == (2O - 1)!!
  #  
  #  We therefore want NRooted(A) * NRooted(O + 1) * NRooted(B)
  #  
  #  O = overlap = bigSplit - smallSplit
  #  bigSplit - overlap = smallSplit = either A or B
  #  n - bigSplit = either B or A
  
  
  # Return:
  NRooted(overlap + 1L) * 
    NRooted(smallSplit) *
    NRooted(n - bigSplit)
}

#' @describeIn SplitMutualInformation Natural logarithm of `TreesConsistentWithTwoSplits`.
#' @importFrom TreeSearch LogTreesMatchingSplit LnRooted.int
#' @export
LogTreesConsistentWithTwoSplits <- function (n, A1, A2=A1) {
  smallSplit <- min(A1, A2)
  bigSplit <- max(A1, A2)
  
  if (smallSplit == 0) return (LogTreesMatchingSplit(bigSplit, n - bigSplit))
  if (bigSplit == n) return (LogTreesMatchingSplit(smallSplit, n - smallSplit))
  
  # Return:
  LnRooted.int(bigSplit - smallSplit + 1L) + 
    LnRooted.int(smallSplit) + 
    LnRooted.int(n - bigSplit)
}
