library("TreeTools")
ReduceTrees <- TreeDist::ReduceTrees

Tree <- function(txt) ape::read.tree(text = txt)
AsTips <- function(sp) {
  which(as.logical(if (TipsInSplits(sp) > nTip / 2) !sp else sp))
}

nTip <- 9

# There are six UnrootedTreeShapes for 9 leaves
# plot(UnrootedTreeWithKey(0, 9))
shape0 <- Tree("((c1, c2), (s, (t, (u, (v, (w, (h1, h2)))))));")
sp0 <- as.Splits(shape0)
tis <- TipsInSplits(sp0)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

fours <- unname(which(tiss == 4))
trios <- unname(which(tiss == 3))
pairs <- seq_len(nTip - 3)[-c(trios, fours)]

mid1 <- xor(sp0[[fours[[1]]]], sp0[[trios[[1]]]])
if (TipsInSplits(mid1) %in% c(1, nTip - 1)) {
  mid2 <- xor(sp0[[fours[[2]]]], sp0[[trios[[2]]]])
} else {
  trios <- trios[2:1]
  mid1 <- xor(sp0[[fours[[1]]]], sp0[[trios[[1]]]])
  mid2 <- xor(sp0[[fours[[2]]]], sp0[[trios[[2]]]])
}

# Align trio1 with mid1
if (!TipsInSplits(xor(sp0[[trios[[1]]]], sp0[[pairs[[1]]]])) %in% 
    c(1, nTip - 1)) {
  pairs <- pairs[2:1]
}

canonOrder <- TipLabels(shape0)[c(
  centre = AsTips(xor(sp0[[fours[[1]]]], sp0[[fours[[2]]]])),
  mid1 = AsTips(mid1),
  trio1 = AsTips(xor(sp0[[trios[[1]]]], sp0[[pairs[[1]]]])),
  pair1 = AsTips(sp0[[pairs[[1]]]]),
  mid2 = AsTips(mid2),
  trio2 = AsTips(xor(sp0[[trios[[2]]]], sp0[[pairs[[2]]]])),
  pair2 = AsTips(sp0[[pairs[[2]]]])
  )]
shape0 <- RenumberTips(shape0, canonOrder)

trees0 <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
scores0 <- sapply(seq_along(trees0), function(i) {
  reduced <- ReduceTrees(shape0, trees0[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
valid0 <- !is.na(scores0)

splits0 <- vapply(which(valid0), function(i) {
  as.integer(!(trees0[[i]] |> as.Splits() |> PolarizeSplits(nTip))) |> sort()
}, integer(nTip - 3))


shape1 <- Tree("((c1, c2), (s, (t, (u, ((p1, p2), (q1, q2))))));")
sp1 <- as.Splits(shape1)
tis <- TipsInSplits(sp1)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

fours <- unname(which(tiss == 4))
trio <- unname(which.max(tiss == 3))
pairs <- seq_len(nTip - 3)[-c(trio, fours)]

trioSp <- sp1[[trio]]
for (trioPair in pairs) {
  soloSp <- xor(trioSp, sp1[[trioPair]])
  if (TipsInSplits(soloSp) %in% c(1, nTip - 1)) break
}

if (!TipsInSplits(xor(trioSp, sp1[[fours[[1]]]])) %in% c(1, nTip - 1)) {
  fours <- fours[2:1]
}

otherPairs <- pairs[pairs != trioPair]

canonOrder <- TipLabels(shape1)[c(
  centre = AsTips(xor(sp1[[fours[[1]]]], sp1[[trio]])),
  mid1 = AsTips(soloSp),
  cherry = AsTips(sp1[[trioPair]]),
  mid2 = AsTips(xor(sp1[[fours[[2]]]], sp1[[fours[[1]]]])),
  pair1 = AsTips(sp1[[otherPairs[[1]]]]),
  pair2 = AsTips(sp1[[otherPairs[[2]]]])
  )]
shape1 <- RenumberTips(shape1, canonOrder)

trees1 <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
scores1 <- sapply(seq_along(trees1), function(i) {
  reduced <- ReduceTrees(shape1, trees1[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
valid1 <- !is.na(scores1)

splits1 <- vapply(which(valid1), function(i) {
  as.integer(!(trees1[[i]] |> as.Splits() |> PolarizeSplits(nTip))) |> sort()
}, integer(nTip - 3))



shape2 <- Tree("((c1, c2), (s, (t, ((p1, p2), (u, (q1, q2))))));")
sp2 <- as.Splits(shape2)
tis <- TipsInSplits(sp2)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

fours <- unname(which.max(tiss == 4))
trios <- unname(which(tiss == 3))
pairs <- seq_len(nTip - 3)[-c(trios, fours)]

if (TipsInSplits(xor(sp2[[trios[[2]]]], sp2[[fours]])) %in% c(1, nTip - 1)) {
  trios <- trios[2:1]
}

trioSp <- sp2[[trios[[1]]]]
for (trioPair1 in pairs) {
  soloSp <- xor(trioSp, sp2[[trioPair1]])
  if (TipsInSplits(soloSp) %in% c(1, nTip - 1)) break
}
trioSp2 <- sp2[[trios[[2]]]]
for (trioPair2 in setdiff(pairs, trioPair1)) {
  soloSp2 <- xor(trioSp2, sp2[[trioPair2]])
  if (TipsInSplits(soloSp2) %in% c(1, nTip - 1)) break
}

canonOrder <- TipLabels(shape2)[c(
  s = AsTips(soloSp),
  c = AsTips(sp2[[trioPair1]]),
  t = AsTips(xor(sp2[[fours]], trioSp)),
  u = AsTips(soloSp2),
  q = AsTips(sp2[[trioPair2]]),
  p = AsTips(sp2[[setdiff(pairs, c(trioPair1, trioPair2))]])
)]
shape2 <- RenumberTips(shape2, canonOrder)

trees2 <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
scores2 <- sapply(seq_along(trees2), function(i) {
  reduced <- ReduceTrees(shape2, trees2[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
valid2 <- !is.na(scores2)

splits2 <- vapply(which(valid2), function(i) {
  as.integer(!(trees2[[i]] |> as.Splits() |> PolarizeSplits(nTip))) |> sort()
}, integer(nTip - 3))


shape3 <- Tree("((c1, c2), (s, ((h1, h2), ((p1, p2), (q1, q2)))));")
sp3 <- as.Splits(shape3)
tis <- TipsInSplits(sp3)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

fours <- unname(which(tiss == 4))
trios <- unname(which(tiss == 3))
pairs <- seq_len(nTip - 3)[-c(trios, fours)]

trioSp <- sp3[[trios]]
for (trioPair in pairs) {
  soloSp <- xor(trioSp, sp3[[trioPair]])
  if (TipsInSplits(soloSp) %in% c(1, nTip - 1)) break
}
soloTip <- AsTips(soloSp)
quadSp <- PolarizeSplits(sp3[[fours]], soloTip)

for (midPair in pairs[pairs != trioPair]) {
  midPairSp <- quadSp & PolarizeSplits(sp3[[midPair]], soloTip)
  if (TipsInSplits(midPairSp) == 3) break
}

pairs <- setdiff(pairs, c(trioPair, midPair))

canonOrder <- TipLabels(shape3)[c(
  solo = soloTip,
  soloPair = AsTips(sp3[[trioPair]]),
  midPair = AsTips(sp3[[midPair]]),
  pair1 = AsTips(sp3[[pairs[[1]]]]),
  pair2 = AsTips(sp3[[pairs[[2]]]])
)]
shape3 <- RenumberTips(shape3, canonOrder)

trees3 <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
scores3 <- sapply(seq_along(trees3), function(i) {
  reduced <- ReduceTrees(shape3, trees3[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
valid3 <- !is.na(scores3)

splits3 <- vapply(which(valid3), function(i) {
  as.integer(!(trees3[[i]] |> as.Splits() |> PolarizeSplits(nTip))) |> sort()
}, integer(nTip - 3))



shape4 <- Tree("((c1, c2), (s, ((t, (p1, p2)), (u, (q1, q2)))));")
sp4 <- as.Splits(shape4)
tis <- TipsInSplits(sp4)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

trios <- unname(which(tiss == 3))
pairs <- seq_len(nTip - 3)[-trios]

trioSp1 <- sp4[[trios[[1]]]]
for (trioPair1 in pairs) {
  soloSp1 <- xor(trioSp1, sp4[[trioPair1]])
  if (TipsInSplits(soloSp1) %in% c(1, nTip - 1)) break
}
trioSp2 <- sp4[[trios[[2]]]]
for (trioPair2 in setdiff(pairs, trioPair1)) {
  soloSp2 <- xor(trioSp2, sp4[[trioPair2]])
  if (TipsInSplits(soloSp2) %in% c(1, nTip - 1)) break
}
trioPair3 <- pairs[!pairs %in% c(trioPair1, trioPair2)]
trioSp3 <- sp4[[trios[[3]]]]
soloSp3 <- xor(trioSp3, sp4[[trioPair3]])

canonOrder <- TipLabels(shape4)[c(
  solo1 = AsTips(soloSp1),
  pair1 = AsTips(sp4[[trioPair1]]),
  solo2 = AsTips(soloSp2),
  pair2 = AsTips(sp4[[trioPair2]]),
  solo3 = AsTips(soloSp3),
  pair3 = AsTips(sp4[[trioPair3]])
)]
shape4 <- RenumberTips(shape4, canonOrder)

trees4 <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
scores4 <- sapply(seq_along(trees4), function(i) {
  reduced <- ReduceTrees(shape4, trees4[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
valid4 <- !is.na(scores4)

splits4 <- vapply(which(valid4), function(i) {
  as.integer(!(trees4[[i]] |> as.Splits() |> PolarizeSplits(nTip))) |> sort()
}, integer(nTip - 3))




shape5 <- Tree("(((c1, c2), (h1, h2)), (s, ((p1, p2)), (q1, q2)));")
sp5 <- as.Splits(shape5)
tis <- TipsInSplits(sp5)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

fours <- unname(which(tiss == 4))
trios <- unname(which(tiss == 3))
pairs <- seq_len(nTip - 3)[-c(trios, fours)]

mid1 <- xor(sp5[[fours[[1]]]], sp5[[trios[[1]]]])
if (TipsInSplits(mid1) %in% c(1, nTip - 1)) {
  mid2 <- xor(sp5[[fours[[2]]]], sp5[[trios[[2]]]])
} else {
  trios <- trios[2:1]
  mid1 <- xor(sp5[[fours[[1]]]], sp5[[trios[[1]]]])
  mid2 <- xor(sp5[[fours[[2]]]], sp5[[trios[[2]]]])
}

# Align trio1 with mid1
if (!TipsInSplits(xor(sp5[[trios[[1]]]], sp5[[pairs[[1]]]])) %in% 
    c(1, nTip - 1)) {
  pairs <- pairs[2:1]
}

canonOrder <- TipLabels(shape5)[c(
  centre = AsTips(xor(sp5[[fours[[1]]]], sp5[[fours[[2]]]])),
  mid1 = AsTips(mid1),
  trio1 = AsTips(xor(sp5[[trios[[1]]]], sp5[[pairs[[1]]]])),
  pair1 = AsTips(sp5[[pairs[[1]]]]),
  mid2 = AsTips(mid2),
  trio2 = AsTips(xor(sp5[[trios[[2]]]], sp5[[pairs[[2]]]])),
  pair2 = AsTips(sp5[[pairs[[2]]]])
)]
shape5 <- RenumberTips(shape5, canonOrder)

trees5 <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
scores5 <- sapply(seq_along(trees5), function(i) {
  reduced <- ReduceTrees(shape5, trees5[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
valid5 <- !is.na(scores5)

splits5 <- vapply(which(valid5), function(i) {
  as.integer(!(trees5[[i]] |> as.Splits() |> PolarizeSplits(nTip))) |> sort()
}, integer(nTip - 3))




# Define packing algorithm based on range
library("bit64")
offset <- c(
  min(splits1[1, ], splits2[1, ], splits3[1, ],
      splits4[1, ], splits5[1, ], splits6[1, ]),
  min(splits1[2, ], splits2[2, ], splits3[2, ],
      splits4[2, ], splits5[2, ], splits6[2, ]),
  min(splits1[3, ], splits2[3, ], splits3[3, ],
      splits4[3, ], splits5[3, ], splits6[3, ]),
  min(splits1[4, ], splits2[4, ], splits3[4, ],
      splits4[4, ], splits5[4, ], splits6[4, ]),
  min(splits1[5, ], splits2[5, ], splits3[5, ],
      splits4[5, ], splits5[5, ], splits6[5, ]),
  min(splits1[6, ], splits2[6, ], splits3[6, ],
      splits4[6, ], splits5[6, ], splits6[6, ])
  ) |>
  as.integer64()

rng <-  c(
  max(splits1[1, ], splits2[1, ], splits3[1, ],
      splits4[1, ], splits5[1, ], splits6[1, ]),
  max(splits1[2, ], splits2[2, ], splits3[2, ],
      splits4[2, ], splits5[2, ], splits6[2, ]),
  max(splits1[3, ], splits2[3, ], splits3[3, ],
      splits4[3, ], splits5[3, ], splits6[3, ]),
  max(splits1[4, ], splits2[4, ], splits3[4, ],
      splits4[4, ], splits5[4, ], splits6[4, ]),
  max(splits1[5, ], splits2[5, ], splits3[5, ],
      splits4[5, ], splits5[5, ], splits6[5, ]),
  max(splits1[6, ], splits2[6, ], splits3[6, ],
      splits4[6, ], splits5[6, ], splits6[6, ])
) -  c(
  min(splits1[1, ], splits2[1, ], splits3[1, ],
      splits4[1, ], splits5[1, ], splits6[1, ]),
  min(splits1[2, ], splits2[2, ], splits3[2, ],
      splits4[2, ], splits5[2, ], splits6[2, ]),
  min(splits1[3, ], splits2[3, ], splits3[3, ],
      splits4[3, ], splits5[3, ], splits6[3, ]),
  min(splits1[4, ], splits2[4, ], splits3[4, ],
      splits4[4, ], splits5[4, ], splits6[4, ]),
  min(splits1[5, ], splits2[5, ], splits3[5, ],
      splits4[5, ], splits5[5, ], splits6[5, ]),
  min(splits1[6, ], splits2[6, ], splits3[6, ],
      splits4[6, ], splits5[6, ], splits6[6, ])
)

BitPack9 <- function(vec) {
  v <- as.integer64(vec)
  as.character(
    (v[1] - offset[[1]]) * 134217728L +
      (v[2] - offset[[2]]) * 1048576L +
      (v[3] - offset[[3]]) * 8192L +
      (v[4] - offset[[4]]) * 64L +
      (v[5] - offset[[5]]) * 64L +
      (v[6] - offset[[6]]))
}

pecPack <- apply(pecSplits, 2, BitPack8)
pecDF <- data.frame(key = pecPack, score = pecScores[pecValid])
pecDF <- pecDF[order(pecDF$key), ]

mixPack <- apply(mixSplits, 2, BitPack8)
mixDF <- data.frame(key = mixPack, score = mixScores[mixValid])
mixDF <- mixDF[order(mixDF$key), ]

midPack <- apply(midSplits, 2, BitPack8)
midDF <- data.frame(key = midPack, score = midScores[midValid])
midDF <- midDF[order(midDF$key), ]

balPack <- apply(balSplits, 2, BitPack8)
balDF <- data.frame(key = balPack, score = balScores[balValid])
balDF <- balDF[order(balDF$key), ]


header_content <- paste0(
  "// Generated from data-raw/spr-exact.R\n",
  "#include <cstdint>\n#include <array>\n#include <algorithm>\n\n",
  "struct SPRScore64 { uint64_t key; int score; };\n\n",
  
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
