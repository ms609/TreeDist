library("TreeTools")
if (!exists("ReduceTrees")) ReduceTrees <- TreeDist::ReduceTrees

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
saveRDS(pecScores, "pecScores.rds")
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

quad <- which(tiss == 4)
trio <- which(tiss == 3)
pairs <- (1:5)[-c(quad, trio)]

trioSp <- mixSp[[trio]]
sTip <- xor(mixSp[[quad]], trioSp)
for (trioPair in pairs) {
  soloSp <- xor(trioSp, mixSp[[trioPair]])
  if (TipsInSplits(soloSp) %in% c(1, nTip - 1)) break
}
cherries <- pairs[pairs != trioPair]

.FewerTips <- function(sp) {
  which(as.logical(if (TipsInSplits(sp) > nTip / 2) !sp else sp))
}
midTip <- .FewerTips(sTip)
trioTip <- .FewerTips(xor(mixSp[[trioPair]], trioSp))
trioPairTip <- .FewerTips(mixSp[[trioPair]])
otherSp1 <- .FewerTips(mixSp[[cherries[[1]]]])
otherSp2 <- .FewerTips(mixSp[[cherries[[2]]]])
canonOrder <- TipLabels(mix)[
  c(midTip, trioTip, trioPairTip, otherSp1, otherSp2)
  ]

mix <- RenumberTips(mix, canonOrder)

mixTrees <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
mixScores <- vapply(seq_along(mixTrees), function(i) {
  reduced <- ReduceTrees(mix, mixTrees[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA_real_)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
}, double(1))
saveRDS(mixScores, "mixScores.rds")

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

trios <- unname(which(tiss == 3))
trioSp1 <- midSp[[trios[[1]]]]
for (trioPair1 in seq_along(tis)[-trios]) {
  solo1 <- xor(midSp[[trioPair1]], trioSp1)
  if (TipsInSplits(solo1) %in% c(1, nTip - 1)) break
}

trioSp2 <- midSp[[trios[[2]]]]
for (trioPair2 in seq_along(tis)[-c(trios, trioPair1)]) {
  solo2 <- xor(midSp[[trioPair2]], trioSp2)
  if (TipsInSplits(solo2) %in% c(1, nTip - 1)) break
}

canonOrder <- TipLabels(mid)[c(
  .FewerTips(solo1),
  .FewerTips(midSp[[trioPair1]]),
  .FewerTips(solo2),
  .FewerTips(midSp[[trioPair2]]),
  .FewerTips(midSp[[setdiff(1:5, c(trios, trioPair1, trioPair2))]])
)]

mid <- RenumberTips(mid, canonOrder)

midTrees <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
midScores <- vapply(seq_along(midTrees), function(i) {
  reduced <- ReduceTrees(mid, midTrees[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA_real_)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
}, double(1))
saveRDS(midScores, "midScores.rds")

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

cherries <- tiss == 2

canonOrder <- TipLabels(bal)[unlist(
  lapply(which(cherries), function(idx) .FewerTips(balSp[[idx]])),
  use.names = FALSE)
]
bal <- RenumberTips(bal, canonOrder)

balTrees <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
balScores <- vapply(seq_along(balTrees), function(i) {
  reduced <- ReduceTrees(bal, balTrees[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA_real_)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
}, double(1))
saveRDS(balScores, "balScores.rds")

balValid <- !is.na(balScores)

balSplits <- vapply(which(balValid), function(i) {
  as.integer(!(balTrees[[i]] |> as.Splits() |> PolarizeSplits(nTip))) |>
    sort()
}, integer(nTip - 3))

# Define packing algorithm based on range
library("bit64")
offset <- c(
  min(pecSplits[1, ], mixSplits[1, ], midSplits[1, ], balSplits[1, ]),
  min(pecSplits[2, ], mixSplits[2, ], midSplits[2, ], balSplits[2, ]),
  min(pecSplits[3, ], mixSplits[3, ], midSplits[3, ], balSplits[3, ]),
  min(pecSplits[4, ], mixSplits[4, ], midSplits[4, ], balSplits[4, ]),
  min(pecSplits[5, ], mixSplits[5, ], midSplits[5, ], balSplits[5, ]))
print(offset)
offset <- as.integer64(offset)

BitPack8 <- function(vec) {
  v <- as.integer64(vec)
  as.character(
    (v[1] - offset[[1]]) * 134217728L +
      (v[2] - offset[[2]]) * 1048576L +
      (v[3] - offset[[3]]) * 8192L +
      (v[4] - offset[[4]]) * 64L +
      (v[5] - offset[[5]]))
}

pecPack <- apply(pecSplits, 2, BitPack8)
pecDF <- data.frame(key = pecPack, score = pecScores[pecValid])
pecDF <- pecDF[order(sprintf("%020s", pecDF$key)), ]

mixPack <- apply(mixSplits, 2, BitPack8)
mixDF <- data.frame(key = mixPack, score = mixScores[mixValid])
mixDF <- mixDF[order(sprintf("%020s", mixDF$key)), ]

midPack <- apply(midSplits, 2, BitPack8)
midDF <- data.frame(key = midPack, score = midScores[midValid])
midDF <- midDF[order(sprintf("%020s", midDF$key)), ]

balPack <- apply(balSplits, 2, BitPack8)
balDF <- data.frame(key = balPack, score = balScores[balValid])
balDF <- balDF[order(sprintf("%020s", balDF$key)), ]


header_content <- paste0(
  "// Generated from data-raw/spr-exact.R\n",
  "#include <cstdint>\n#include <array>\n#include <algorithm>\n",
  "#include \"lookup.h\"\n\n",
  
  "static constexpr std::array<SPRScore64, ", nrow(pecDF), "> PEC_LOOKUP",
  nTip, " = {{\n",
  paste0("    {", pecDF$key, "ULL, ", pecDF$score, "}", collapse = ",\n"),
  "\n}};\n",
  
  "static constexpr std::array<SPRScore64, ", nrow(mixDF), "> MIX_LOOKUP",
  nTip, " = {{\n",
  paste0("    {", mixDF$key, "ULL, ", mixDF$score, "}", collapse = ",\n"),
  "\n}};\n",
  
  "static constexpr std::array<SPRScore64, ", nrow(midDF), "> MID_LOOKUP",
  nTip, " = {{\n",
  paste0("    {", midDF$key, "ULL, ", midDF$score, "}", collapse = ",\n"),
  "\n}};\n",
  
  "static constexpr std::array<SPRScore64, ", nrow(balDF), "> BAL_LOOKUP",
  nTip, " = {{\n",
  paste0("    {", balDF$key, "ULL, ", balDF$score, "}", collapse = ",\n"),
  "\n}};"
)

writeLines(header_content, sprintf("src/spr/lookup_table_%d.h", nTip))
