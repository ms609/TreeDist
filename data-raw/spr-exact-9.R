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
if (file.exists("scores0.rds")) {
  scores0 <- readRDS("scores0.rds")
} else {
  scores0 <- sapply(cli::cli_progress_along(trees0), function(i) {
    reduced <- ReduceTrees(shape0, trees0[[i]])
    r1 <- reduced[[1]]
    if (is.null(r1) || NTip(r1) != nTip) return(NA)
    r2 <- reduced[[2]]
    TBRDist::USPRDist(r1, r2)
  })
  saveRDS(scores0, file = "scores0.rds")
}
valid0 <- !is.na(scores0)
  
if (file.exists("splits0.rds")) {
  splits0 <- readRDS("splits0.rds")
} else {
  splits0 <- vapply(which(valid0), function(i) {
    as.integer(!(trees0[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
  }, integer(nTip - 3))
  saveRDS(splits0, file = "splits0.rds")
}


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
if (file.exists("scores1.rds")) {
  scores1 <- readRDS("scores1.rds")
} else {
  scores1 <- sapply(cli::cli_progress_along(trees1), function(i) {
    reduced <- ReduceTrees(shape1, trees1[[i]])
    r1 <- reduced[[1]]
    if (is.null(r1) || NTip(r1) != nTip) return(NA)
    r2 <- reduced[[2]]
    TBRDist::USPRDist(r1, r2)
  })
  saveRDS(scores1, file = "scores1.rds")
}
valid1 <- !is.na(scores1)
  
if (file.exists("splits1.rds")) {
  splits1 <- readRDS("splits1.rds")
} else {
  splits1 <- vapply(which(valid1), function(i) {
    as.integer(!(trees1[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
  }, integer(nTip - 3))
  saveRDS(splits1, file = "splits1.rds")
}


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
if (file.exists("scores2.rds")) {
  scores2 <- readRDS("scores2.rds")
} else {
  scores2 <- sapply(cli::cli_progress_along(trees2), function(i) {
    reduced <- ReduceTrees(shape2, trees2[[i]])
    r1 <- reduced[[1]]
    if (is.null(r1) || NTip(r1) != nTip) return(NA)
    r2 <- reduced[[2]]
    TBRDist::USPRDist(r1, r2)
  })
  saveRDS(scores2, file = "scores2.rds")
}
valid2 <- !is.na(scores2)
  
if (file.exists("splits2.rds")) {
  splits2 <- readRDS("splits2.rds")
} else {
  splits2 <- vapply(which(valid2), function(i) {
    as.integer(!(trees2[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
  }, integer(nTip - 3))
  saveRDS(splits2, file = "splits2.rds")
}


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
if (file.exists("scores3.rds")) {
  scores3 <- readRDS("scores3.rds")
} else {
  scores3 <- sapply(cli::cli_progress_along(trees3), function(i) {
    reduced <- ReduceTrees(shape3, trees3[[i]])
    r1 <- reduced[[1]]
    if (is.null(r1) || NTip(r1) != nTip) return(NA)
    r2 <- reduced[[2]]
    TBRDist::USPRDist(r1, r2)
  })
  saveRDS(scores3, file = "scores3.rds")
}
valid3 <- !is.na(scores3)
  
if (file.exists("splits3.rds")) {
  splits3 <- readRDS("splits3.rds")
} else {
  splits3 <- vapply(which(valid3), function(i) {
    as.integer(!(trees3[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
  }, integer(nTip - 3))
  saveRDS(splits3, file = "splits3.rds")
}




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
if (file.exists("scores4.rds")) {
  scores4 <- readRDS("scores4.rds")
} else {
    scores4 <- sapply(cli::cli_progress_along(trees4), function(i) {
    reduced <- ReduceTrees(shape4, trees4[[i]])
    r1 <- reduced[[1]]
    if (is.null(r1) || NTip(r1) != nTip) return(NA)
    r2 <- reduced[[2]]
    TBRDist::USPRDist(r1, r2)
  })
  saveRDS(scores4, file = "scores4.rds")
}
valid4 <- !is.na(scores4)

if (file.exists("splits4.rds")) {
  splits4 <- readRDS("splits4.rds")
} else {
  splits4 <- vapply(which(valid4), function(i) {
    as.integer(!(trees4[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
  }, integer(nTip - 3))
  saveRDS(splits4, file = "splits4.rds")
}




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
if (file.exists("scores5.rds")) {
  scores5 <- readRDS("scores5.rds")
} else {
  scores5 <- sapply(cli::cli_progress_along(trees5), function(i) {
    reduced <- ReduceTrees(shape5, trees5[[i]])
    r1 <- reduced[[1]]
    if (is.null(r1) || NTip(r1) != nTip) return(NA)
    r2 <- reduced[[2]]
    TBRDist::USPRDist(r1, r2)
  })
  saveRDS(scores5, file = "scores5.rds")
}
valid5 <- !is.na(scores5)

if (file.exists("splits5.rds")) {
  splits5 <- readRDS("splits5.rds")
} else {
  splits5 <- vapply(which(valid5), function(i) {
    as.integer(!(trees5[[i]] |> as.Splits() |> PolarizeSplits(nTip))[, 1]) |> sort()
  }, integer(nTip - 3))
  saveRDS(splits5, file = "splits5.rds")
}



source("data-raw/spr-exact-builders.R")

header_content <- paste0(
  c("// Generated in data-raw/spr-exact-9.R",
    DecisionTreeLine("S9_0", scores0, valid0, splits0),
    DecisionTreeLine("S9_1", scores1, valid1, splits1),
    DecisionTreeLine("S9_2", scores2, valid2, splits2),
    DecisionTreeLine("S9_3", scores3, valid3, splits3),
    DecisionTreeLine("S9_4", scores4, valid4, splits4),
    DecisionTreeLine("S9_5", scores5, valid5, splits5)
  )
)

writeLines(header_content, "src/spr/lookup9.h")
