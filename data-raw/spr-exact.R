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
pecDF <- data.frame(key = pecPack, score = pecScores[pecValid])
pecDF <- pecDF[order(pecDF$key), ]

balPack <- apply(balSplits, 2, BitPack7)
balDF <- data.frame(key = balPack, score = balScores[balValid])
balDF <- balDF[order(balDF$key), ]


header_content <- paste0(
  "// Generated from data-raw/spr-exact.R\n",
  "#include <cstdint>\n#include <array>\n#include <algorithm>\n\n",
  "struct SPRScore { uint32_t key; int score; };\n\n",
  "static constexpr std::array<SPRScore, ", nrow(pecDF), "> PEC_LOOKUP = {{\n",
  paste0("    {", pecDF$key, "u, ", pecDF$score, "}", collapse = ",\n"),
  "\n}};\n",
  "static constexpr std::array<SPRScore, ", nrow(balDF), "> BAL_LOOKUP = {{\n",
  paste0("    {", balDF$key, "u, ", balDF$score, "}", collapse = ",\n"),
  "\n}};"
)

writeLines(header_content, "src/spr/lookup_table_7.h")
