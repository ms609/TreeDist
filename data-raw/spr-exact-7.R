library("TreeTools")

Tree <- function(txt) ape::read.tree(text = txt)

pec7 <- Tree("((c1, c2), (s, (t, (u, (h1, h2)))));")
pecSp <- as.Splits(pec7)

nTip <- 7

tis <- TipsInSplits(pecSp)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

trios <- which(tiss == 3)
pairs <- (1:4)[-trios]

midpoint <- xor(pecSp[[trios[[1]]]], pecSp[[trios[[2]]]])
# Align pair1 with trio1
if (!TipsInSplits(xor(pecSp[[pairs[[1]]]], pecSp[[trios[[1]]]])) %in% c(1, 6)) {
  pairs <- pairs[2:1]
}

midTip <- which(as.logical(midpoint))
trio1Tip <- as.logical(xor(pecSp[[trios[[1]]]], pecSp[[pairs[[1]]]]))
if (sum(trio1Tip) == 6) trio1Tip <- !trio1Tip
trio2Tip <- as.logical(xor(pecSp[[trios[[2]]]], pecSp[[pairs[[2]]]]))
if (sum(trio2Tip) == 6) trio2Tip <- !trio2Tip

duo1Tips <- as.logical(pecSp[[pairs[[1]]]])
if (sum(duo1Tips) == 5) duo1Tips <- !duo1Tips
duo2Tips <- as.logical(pecSp[[pairs[[2]]]])
if (sum(duo2Tips) == 5) duo2Tips <- !duo2Tips

canonOrder <- TipLabels(pec7)[
  c(midTip, which(trio1Tip), which(duo1Tips), which(trio2Tip), which(duo2Tips))
]
pec7 <- RenumberTips(pec7, canonOrder)

pecTrees <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
pecScores <- sapply(seq_along(pecTrees), function(i) {
  reduced <- ReduceTrees(pec7, pecTrees[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
pecValid <- !is.na(pecScores)

pecSplits <- vapply(which(pecValid), function(i) {
  as.integer(!(pecTrees[[i]] |> as.Splits() |> PolarizeSplits(7))) |> sort()
}, integer(4))


bal7 <- Tree("(((p1, p2), (q1, q2)), (s, (r1, r2)));")
balSp <- as.Splits(bal7)
tis <- TipsInSplits(balSp)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

firstTrio <- which.max(tiss)
firstTrioSp <- balSp[[firstTrio]]
for (trioPair in seq_along(tis)[-firstTrio]) {
  soloSp <- xor(balSp[[trioPair]], firstTrioSp)
  if (TipsInSplits(soloSp) == nTip - 1) break
}
otherSp <- seq_along(tis)[-c(trioPair, firstTrio)]
singleton <- !soloSp
singleTip <- which(as.logical(singleton))
trioPairTip <- as.logical(balSp[[trioPair]])
if (tisBig[[trioPair]]) trioPairTip <- !as.logical(trioPairTip)
otherSp1 <- as.logical(balSp[[otherSp[[1]]]])
if (tisBig[[otherSp[[1]]]]) otherSp1 <- !as.logical(otherSp1)
otherSp2 <- as.logical(balSp[[otherSp[[2]]]])
if (tisBig[[otherSp[[2]]]]) otherSp2 <- !as.logical(otherSp2)
canonOrder <- TipLabels(bal7)[
  c(singleTip, which(trioPairTip), which(otherSp1), which(otherSp2))]

bal7 <- RenumberTips(bal7, canonOrder)

balTrees <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
balScores <- vapply(seq_along(balTrees), function(i) {
  reduced <- ReduceTrees(bal7, balTrees[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA_real_)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
}, double(1))

balValid <- !is.na(balScores)

balSplits <- vapply(which(balValid), function(i) {
  as.integer(!(balTrees[[i]] |> as.Splits() |> PolarizeSplits(7))) |>
    sort()
}, integer(4))

# Define packing algorithm based on range
range(pecSplits[1, ], balSplits[1, ])
range(pecSplits[2, ], balSplits[2, ])
range(pecSplits[3, ], balSplits[3, ])
range(pecSplits[4, ], balSplits[4, ])
BitPack7 <- function(vec) {
  bitwShiftL(vec[1] - 3, 18) +
    bitwShiftL(vec[2] - 7, 12) +
    bitwShiftL(vec[3] - 15, 6)  +
    vec[4] - 33
}

pecPack <- apply(pecSplits, 2, BitPack7)

.EntropyMap <- function(scores, tags) {
  if (length(scores) == 0) return(NULL)
  h0 <- Ntropy(table(scores))
  if (h0 == 0) return(scores[[1]])
  hJ <- apply(tags, 2, function(sp) {
    Ntropy(table(scores, sp))
  })
  flag <- which.min(hJ)
  if (hJ[flag] == h0) {
    ayeScore <- scores[flag]
    nayScore <- setdiff(scores, ayeScore)
    list(
      q = names(flag),
      aye = ayeScore,
      nay = nayScore
    )
  } else {
    has <- tags[, flag]
    ayeTags <- tags[has, -flag, drop = FALSE]
    nayTags <- tags[!has, -flag, drop = FALSE]
    list(
      q = names(flag),
      aye = .EntropyMap(scores[has], ayeTags[, colSums(ayeTags) > 0, drop = FALSE]),
      nay = .EntropyMap(scores[!has], nayTags[, colSums(nayTags) > 0, drop = FALSE])
    )
  }
}
.Translate <- function(node) {
  if (is.list(node)) {
    sprintf("f(%s,%s,%s)", sub("sp", "", node[[1]]),
            .Translate(node[[2]]),
            .Translate(node[[3]]))
  } else {
    as.character(node)
  }
}

balDF <- data.frame(key = balPack, score = balScores[balValid])
balDF <- balDF[order(balDF$key), ]

KeyEntry <- function(str, df) {
  paste0("alignas(64) static constexpr std::array<uint32_t, ", nrow(df), "> ",
         str, "_KEY7 = {", paste(df$key, collapse = "U,"), "U};")
}
ValEntry <- function(str, df) {
  paste0("alignas(64) static constexpr std::array<uint8_t, ", nrow(df), "> ",
         str, "_VAL7 = {", paste(df$score, collapse = ","), "};")
}
Entries <- function(str, df) {
  c(KeyEntry(str, df), ValEntry(str, df))
}
  
header_content <- paste0(
  c("// Generated in data-raw/spr-exact-7.R",
  "#include <cstdint>\n#include <array>",
  Entries("PEC", pecDF),
  Entries("BAL", balDF))
)

writeLines(header_content, "src/spr/lookup7.h")
