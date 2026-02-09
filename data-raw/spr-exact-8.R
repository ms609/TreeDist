library("TreeTools")

Tree <- function(txt) ape::read.tree(text = txt)
nTip <- 8

# There are four UnrootedTreeShapes for 8 leaves
pec <- Tree("((c1, c2), (s, (t, (u, (v, (h1, h2))))));")
pecSp <- as.Splits(pec)


tis <- TipsInSplits(pecSp)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

fours <- unname(which(tiss == 4))
trios <- unname(which(tiss == 3))
pairs <- seq_len(nTip - 3)[-c(trios, fours)]

mid1 <- xor(pecSp[[fours]], pecSp[[trios[[1]]]])
mid2 <- xor(pecSp[[fours]], pecSp[[trios[[2]]]])
mid1Tip <- which(as.logical(mid1))
mid2Tip <- which(as.logical(mid2))

# Align trio1 with mid1
if (!TipsInSplits(xor(pecSp[[trios[[1]]]], pecSp[[pairs[[1]]]])) %in% 
    c(1, nTip - 1)) {
  pairs <- pairs[2:1]
}

trio1Tip <- as.logical(xor(pecSp[[trios[[1]]]], pecSp[[pairs[[1]]]]))
if (sum(trio1Tip) == nTip - 1) trio1Tip <- !trio1Tip
trio2Tip <- as.logical(xor(pecSp[[trios[[2]]]], pecSp[[pairs[[2]]]]))
if (sum(trio2Tip) == nTip - 1) trio2Tip <- !trio2Tip

duo1Tips <- as.logical(pecSp[[pairs[[1]]]])
if (sum(duo1Tips) == nTip - 2) duo1Tips <- !duo1Tips
duo2Tips <- as.logical(pecSp[[pairs[[2]]]])
if (sum(duo2Tips) == nTip - 2) duo2Tips <- !duo2Tips

canonOrder <- TipLabels(pec)[
  c(mid1Tip, which(trio1Tip), which(duo1Tips),
    mid2Tip, which(trio2Tip), which(duo2Tips))
]
pec <- RenumberTips(pec, canonOrder)

pecTrees <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
pecScores <- sapply(seq_along(pecTrees), function(i) {
  reduced <- ReduceTrees(pec, pecTrees[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
pecValid <- !is.na(pecScores)

pecSplits <- vapply(which(pecValid), function(i) {
  as.integer(!(pecTrees[[i]] |> as.Splits() |> PolarizeSplits(nTip))) |> sort()
}, integer(nTip - 3))




mix <- Tree("(((p1, p2), (q1, q2)), (s, (t, (c1, c2))));")
mixSp <- as.Splits(mix)
tis <- TipsInSplits(mixSp)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

firstTrio <- which.max(tiss)
firstTrioSp <- mixSp[[firstTrio]]
for (trioPair in seq_along(tis)[-firstTrio]) {
  soloSp <- xor(mixSp[[trioPair]], firstTrioSp)
  if (TipsInSplits(soloSp) == nTip - 1) break
}
otherSp <- seq_along(tis)[-c(trioPair, firstTrio)]
singleton <- !soloSp
singleTip <- which(as.logical(singleton))
trioPairTip <- as.logical(mixSp[[trioPair]])
if (tisBig[[trioPair]]) trioPairTip <- !as.logical(trioPairTip)
otherSp1 <- as.logical(mixSp[[otherSp[[1]]]])
if (tisBig[[otherSp[[1]]]]) otherSp1 <- !as.logical(otherSp1)
otherSp2 <- as.logical(mixSp[[otherSp[[2]]]])
if (tisBig[[otherSp[[2]]]]) otherSp2 <- !as.logical(otherSp2)
canonOrder <- TipLabels(mix)[
  c(singleTip, which(trioPairTip), which(otherSp1), which(otherSp2))]

mix <- RenumberTips(mix, canonOrder)

mixTrees <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
mixScores <- vapply(seq_along(mixTrees), function(i) {
  reduced <- ReduceTrees(mix, mixTrees[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA_real_)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
}, double(1))

mixValid <- !is.na(mixScores)

mixSplits <- vapply(which(mixValid), function(i) {
  as.integer(!(mixTrees[[i]] |> as.Splits() |> PolarizeSplits(nTip))) |>
    sort()
}, integer(nTip - 3))




mid <- Tree("((((p1, p2), p0), (q0, (q1, q2))), (c1, c2));")
midSp <- as.Splits(mid)
tis <- TipsInSplits(midSp)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

firstTrio <- which.max(tiss)
firstTrioSp <- midSp[[firstTrio]]
for (trioPair in seq_along(tis)[-firstTrio]) {
  soloSp <- xor(midSp[[trioPair]], firstTrioSp)
  if (TipsInSplits(soloSp) == nTip - 1) break
}
otherSp <- seq_along(tis)[-c(trioPair, firstTrio)]
singleton <- !soloSp
singleTip <- which(as.logical(singleton))
trioPairTip <- as.logical(midSp[[trioPair]])
if (tisBig[[trioPair]]) trioPairTip <- !as.logical(trioPairTip)
otherSp1 <- as.logical(midSp[[otherSp[[1]]]])
if (tisBig[[otherSp[[1]]]]) otherSp1 <- !as.logical(otherSp1)
otherSp2 <- as.logical(midSp[[otherSp[[2]]]])
if (tisBig[[otherSp[[2]]]]) otherSp2 <- !as.logical(otherSp2)
canonOrder <- TipLabels(mid)[
  c(singleTip, which(trioPairTip), which(otherSp1), which(otherSp2))]

mid <- RenumberTips(mid, canonOrder)

midTrees <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
midScores <- vapply(seq_along(midTrees), function(i) {
  reduced <- ReduceTrees(mid, midTrees[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA_real_)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
}, double(1))

midValid <- !is.na(midScores)

midSplits <- vapply(which(midValid), function(i) {
  as.integer(!(midTrees[[i]] |> as.Splits() |> PolarizeSplits(nTip))) |>
    sort()
}, integer(nTip - 3))




bal <- Tree("(((p1, p2), (q1, q2)), ((s1, s2), (r1, r2)));")
balSp <- as.Splits(bal)
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
canonOrder <- TipLabels(bal)[
  c(singleTip, which(trioPairTip), which(otherSp1), which(otherSp2))]

bal <- RenumberTips(bal, canonOrder)

balTrees <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
balScores <- vapply(seq_along(balTrees), function(i) {
  reduced <- ReduceTrees(bal, balTrees[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA_real_)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
}, double(1))

balValid <- !is.na(balScores)

balSplits <- vapply(which(balValid), function(i) {
  as.integer(!(balTrees[[i]] |> as.Splits() |> PolarizeSplits(nTip))) |>
    sort()
}, integer(nTip - 3))

# Define packing algorithm based on range
range(pecSplits[1, ], mixSplits[1, ], midSplits[1, ], balSplits[1, ])
range(pecSplits[2, ], mixSplits[2, ], midSplits[2, ], balSplits[2, ])
range(pecSplits[3, ], mixSplits[3, ], midSplits[3, ], balSplits[3, ])
range(pecSplits[4, ], mixSplits[4, ], midSplits[4, ], balSplits[4, ])
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
  
  "static constexpr std::array<SPRScore, ", nrow(pecDF), "> PEC_LOOKUP",
  nTip, " = {{\n",
  paste0("    {", pecDF$key, "u, ", pecDF$score, "}", collapse = ",\n"),
  "\n}};\n",
  
  "static constexpr std::array<SPRScore, ", nrow(mixDF), "> MIX_LOOKUP",
  nTip, " = {{\n",
  paste0("    {", mixDF$key, "u, ", mixDF$score, "}", collapse = ",\n"),
  "\n}};\n",
  
  "static constexpr std::array<SPRScore, ", nrow(midDF), "> MID_LOOKUP",
  nTip, " = {{\n",
  paste0("    {", midDF$key, "u, ", midDF$score, "}", collapse = ",\n"),
  "\n}};\n",
  
  "static constexpr std::array<SPRScore, ", nrow(balDF), "> BAL_LOOKUP",
  nTip, " = {{\n",
  paste0("    {", balDF$key, "u, ", balDF$score, "}", collapse = ",\n"),
  "\n}};"
)

writeLines(header_content, sprintf("src/spr/lookup_table_%d.h", nTip))
