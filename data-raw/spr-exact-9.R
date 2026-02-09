library("TreeTools")

Tree <- function(txt) ape::read.tree(text = txt)
FewerTips <- function(sp) {
  which(as.logical(if (TipsInSplits(sp) > nTip / 2) !sp else sp))
}

nTip <- 9

# There are six UnrootedTreeShapes for 9 leaves
# plot(UnrootedTreeWithKey(0, 9))
shape0 <- Tree("((c1, c2), (s, (t, (u, (v, (w, (h1, h2)))))));")
shape1 <- Tree("((c1, c2), (s, (t, (u, ((p1, p2), (q1, q2))))));")
shape2 <- Tree("((c1, c2), (s, (t, ((p1, p2), (u, (q1, q2))))));")
shape3 <- Tree("((c1, c2), (s, ((p1, p2), (t, (u, (q1, q2))))));")
shape4 <- Tree("((c1, c2), (s, ((r1, r1), ((p1, p2), (q1, q2)))));")
shape5 <- Tree("((c1, c2), (s, ((t, (p1, p2)), (u, (q1, q2)))));")


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

