library("TreeTools")
if (!exists("ReduceTrees")) ReduceTrees <- TreeDist::ReduceTrees

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
scores0 <- sapply(cli::cli_progress_along(trees0), function(i) {
  reduced <- ReduceTrees(shape0, trees0[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
saveRDS(scores0, file = "scores0.rds")
valid0 <- !is.na(scores0)

splits0 <- vapply(which(valid0), function(i) {
  as.integer(!(trees0[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
}, integer(nTip - 3))
saveRDS(splits0, file = "splits0.rds")


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
scores1 <- sapply(cli::cli_progress_along(trees1), function(i) {
  reduced <- ReduceTrees(shape1, trees1[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
valid1 <- !is.na(scores1)

splits1 <- vapply(which(valid1), function(i) {
  as.integer(!(trees1[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
}, integer(nTip - 3))

saveRDS(scores1, file = "scores1.rds")
saveRDS(splits1, file = "splits1.rds")


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
scores2 <- sapply(cli::cli_progress_along(trees2), function(i) {
  reduced <- ReduceTrees(shape2, trees2[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
valid2 <- !is.na(scores2)

splits2 <- vapply(which(valid2), function(i) {
  as.integer(!(trees2[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
}, integer(nTip - 3))

saveRDS(scores2, file = "scores2.rds")
saveRDS(splits2, file = "splits2.rds")


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
scores3 <- sapply(cli::cli_progress_along(trees3), function(i) {
  reduced <- ReduceTrees(shape3, trees3[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
valid3 <- !is.na(scores3)

splits3 <- vapply(which(valid3), function(i) {
  as.integer(!(trees3[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
}, integer(nTip - 3))

saveRDS(scores3, file = "scores3.rds")
saveRDS(splits3, file = "splits3.rds")



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
scores4 <- sapply(cli::cli_progress_along(trees4), function(i) {
  reduced <- ReduceTrees(shape4, trees4[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
valid4 <- !is.na(scores4)

splits4 <- vapply(which(valid4), function(i) {
  as.integer(!(trees4[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
}, integer(nTip - 3))

saveRDS(scores4, file = "scores4.rds")
saveRDS(splits4, file = "splits4.rds")




shape5 <- Tree("(((c1, c2), (h1, h2)), (s, ((p1, p2), (q1, q2))));")
sp5 <- as.Splits(shape5)
tis <- TipsInSplits(sp5)
tisBig <- tis > nTip / 2
tiss <- tis
tiss[tisBig] <- nTip - tis[tisBig]

fours <- unname(which(tiss == 4))
pairs <- seq_len(nTip - 3)[-fours]
soloTip <- AsTips(xor(sp5[[fours[[1]]]], sp5[[fours[[2]]]]))
sp5 <- !PolarizeSplits(sp5, soloTip)

four1 <- sp5[[fours[[1]]]]
pairedWith1 <- vapply(pairs, function(p) {
  TipsInSplits(four1 & sp5[[p]]) > 0
}, logical(1))

pairs1 <- pairs[pairedWith1]
pairs2 <- pairs[!pairedWith1]

canonOrder <- TipLabels(shape5)[c(
  solo = soloTip,
  pair1a = AsTips(sp5[[pairs1[[1]]]]),
  pair1b = AsTips(sp5[[pairs1[[2]]]]),
  pair2a = AsTips(sp5[[pairs2[[1]]]]),
  pair2b = AsTips(sp5[[pairs2[[2]]]])
)]
shape5 <- RenumberTips(shape5, canonOrder)

trees5 <- as.phylo(seq_len(NUnrooted(nTip)), nTip, canonOrder)
scores5 <- sapply(cli::cli_progress_along(trees5), function(i) {
  reduced <- ReduceTrees(shape5, trees5[[i]])
  r1 <- reduced[[1]]
  if (is.null(r1) || NTip(r1) != nTip) return(NA)
  r2 <- reduced[[2]]
  TBRDist::USPRDist(r1, r2)
})
valid5 <- !is.na(scores5)

splits5 <- vapply(which(valid5), function(i) {
  as.integer(!(trees5[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
}, integer(nTip - 3))

saveRDS(scores5, file = "scores5.rds")
saveRDS(splits5, file = "splits5.rds")




# Define packing algorithm based on range
library("bit64")
scores0 <- readRDS("scores0.rds")
sp0 <- readRDS("splits0.rds")
sp1 <- splits1
sp2 <- splits2
sp3 <- splits3
sp4 <- splits4
sp5 <- splits5
offset <- c(
  min(sp0[1, ], sp1[1, ], sp2[1, ], sp3[1, ], sp4[1, ], sp5[1, ]),
  min(sp0[2, ], sp1[2, ], sp2[2, ], sp3[2, ], sp4[2, ], sp5[2, ]),
  min(sp0[3, ], sp1[3, ], sp2[3, ], sp3[3, ], sp4[3, ], sp5[3, ]),
  min(sp0[4, ], sp1[4, ], sp2[4, ], sp3[4, ], sp4[4, ], sp5[4, ]),
  min(sp0[5, ], sp1[5, ], sp2[5, ], sp3[5, ], sp4[5, ], sp5[5, ]),
  min(sp0[6, ], sp1[6, ], sp2[6, ], sp3[6, ], sp4[6, ], sp5[6, ])
  ) |>
  as.integer64()

rng <-  c(
  max(sp0[1, ], sp1[1, ], sp2[1, ], sp3[1, ], sp4[1, ], sp5[1, ]),
  max(sp0[2, ], sp1[2, ], sp2[2, ], sp3[2, ], sp4[2, ], sp5[2, ]),
  max(sp0[3, ], sp1[3, ], sp2[3, ], sp3[3, ], sp4[3, ], sp5[3, ]),
  max(sp0[4, ], sp1[4, ], sp2[4, ], sp3[4, ], sp4[4, ], sp5[4, ]),
  max(sp0[5, ], sp1[5, ], sp2[5, ], sp3[5, ], sp4[5, ], sp5[5, ]),
  max(sp0[6, ], sp1[6, ], sp2[6, ], sp3[6, ], sp4[6, ], sp5[6, ])
) -  c(
  min(sp0[1, ], sp1[1, ], sp2[1, ], sp3[1, ], sp4[1, ], sp5[1, ]),
  min(sp0[2, ], sp1[2, ], sp2[2, ], sp3[2, ], sp4[2, ], sp5[2, ]),
  min(sp0[3, ], sp1[3, ], sp2[3, ], sp3[3, ], sp4[3, ], sp5[3, ]),
  min(sp0[4, ], sp1[4, ], sp2[4, ], sp3[4, ], sp4[4, ], sp5[4, ]),
  min(sp0[5, ], sp1[5, ], sp2[5, ], sp3[5, ], sp4[5, ], sp5[5, ]),
  min(sp0[6, ], sp1[6, ], sp2[6, ], sp3[6, ], sp4[6, ], sp5[6, ])
)

2^as.integer64(cumsum(rev(ceiling(log2(rng)))))

BitPack9 <- function(vec) {
  v <- as.integer64(vec)
  as.character(
    (v[1] - offset[[1]]) * as.integer64("549755813888") +
      (v[2] - offset[[2]]) * as.integer64("2147483648") +
      (v[3] - offset[[3]]) * 8388608L +
      (v[4] - offset[[4]]) * 32768L +
      (v[5] - offset[[5]]) * 128L +
      (v[6] - offset[[6]]))
}
.Order <- function(keys) {
  order(sprintf("%020s", keys))
}

pack0 <- apply(sp0, 2, BitPack9)
df0 <- data.frame(key = pack0, score = scores0[valid0])
df0 <- df0[.Order(df0$key), ]

pack1 <- apply(sp1, 2, BitPack9)
df1 <- data.frame(key = pack1, score = scores1[valid1])
df1 <- df1[.Order(df1$key), ]

pack2 <- apply(sp2, 2, BitPack9)
df2 <- data.frame(key = pack2, score = scores2[valid2])
df2 <- df2[.Order(df2$key), ]

pack3 <- apply(sp3, 2, BitPack9)
df3 <- data.frame(key = pack3, score = scores3[valid3])
df3 <- df3[.Order(df3$key), ]

pack4 <- apply(sp4, 2, BitPack9)
df4 <- data.frame(key = pack4, score = scores4[valid4])
df4 <- df4[.Order(df4$key), ]

pack5 <- apply(sp5, 2, BitPack9)
df5 <- data.frame(key = pack5, score = scores5[valid5])
df5 <- df5[.Order(df5$key), ]

header_content <- paste0(
  "// Generated from data-raw/spr-exact.R\n",
  "#include <cstdint>\n#include <array>\n#include <algorithm>\n\n",
  "struct SPRScore64 { uint64_t key; int score; };\n\n",
  
  "static constexpr std::array<SPRScore64, ", nrow(df0), "> LOOKUP9_0",
  " = {{\n",
  paste0("    {", df0$key, "ULL, ", df0$score, "}", collapse = ",\n"),
  "\n}};\n",
  "static constexpr std::array<SPRScore64, ", nrow(df1), "> LOOKUP9_1",
  " = {{\n",
  paste0("    {", df1$key, "ULL, ", df1$score, "}", collapse = ",\n"),
  "\n}};\n",
  "static constexpr std::array<SPRScore64, ", nrow(df2), "> LOOKUP9_2",
  " = {{\n",
  paste0("    {", df2$key, "ULL, ", df2$score, "}", collapse = ",\n"),
  "\n}};\n",
  "static constexpr std::array<SPRScore64, ", nrow(df3), "> LOOKUP9_3",
  " = {{\n",
  paste0("    {", df3$key, "ULL, ", df3$score, "}", collapse = ",\n"),
  "\n}};\n",
  "static constexpr std::array<SPRScore64, ", nrow(df4), "> LOOKUP9_4",
  " = {{\n",
  paste0("    {", df4$key, "ULL, ", df4$score, "}", collapse = ",\n"),
  "\n}};\n",
  "static constexpr std::array<SPRScore64, ", nrow(df5), "> LOOKUP9_5",
  " = {{\n",
  paste0("    {", df5$key, "ULL, ", df5$score, "}", collapse = ",\n"),
  "\n}};\n"
)

writeLines(header_content, sprintf("src/spr/lookup_table_%d.h", nTip))
