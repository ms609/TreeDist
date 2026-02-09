library("TreeTools")
ReduceTrees <- TreeDist::ReduceTrees

Tree <- function(txt) ape::read.tree(text = txt)
FewerTips <- function(sp) {
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
  centre = FewerTips(xor(sp0[[fours[[1]]]], sp0[[fours[[2]]]])),
  mid1 = FewerTips(mid1),
  trio1 = FewerTips(xor(sp0[[trios[[1]]]], sp0[[pairs[[1]]]])),
  pair1 = FewerTips(sp0[[pairs[[1]]]]),
  mid2 = FewerTips(mid2),
  trio2 = FewerTips(xor(sp0[[trios[[2]]]], sp0[[pairs[[2]]]])),
  pair2 = FewerTips(sp0[[pairs[[2]]]])
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
  centre = FewerTips(xor(sp1[[fours[[1]]]], sp1[[trio]])),
  mid1 = FewerTips(soloSp),
  cherry = FewerTips(sp1[[trioPair]]),
  mid2 = FewerTips(xor(sp1[[fours[[2]]]], sp1[[fours[[1]]]])),
  pair1 = FewerTips(sp1[[otherPairs[[1]]]]),
  pair2 = FewerTips(sp1[[otherPairs[[2]]]])
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

