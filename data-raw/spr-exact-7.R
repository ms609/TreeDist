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







.EntropyMap <- function(scores, tags) {
  if (length(unique(scores)) == 1) return(scores[[1]])
  
  hJ <- apply(tags, 2, function(sp) Ntropy(table(scores, sp)))
  flag <- which.min(hJ)
  
  has <- tags[, flag]
  if (all(has)) stop("Splitting on allT - fail")
  if (all(!has)) stop("Splitting on allF - fail")
  
  .Retain <- function(rows) {
    ret <- tags[rows, -flag, drop = FALSE]
    cs <- colSums(ret)
    ret[, cs > 0 & cs < dim(ret)[[1]], drop = FALSE]
  }
  
  list(
    q = names(flag),
    aye = .EntropyMap(scores[has], .Retain(has)),
    nay = .EntropyMap(scores[!has], .Retain(!has))
  )
}

.FlattenMap <- function(node) {
  # 1. First pass: Linearize the tree into a list of nodes
  nodes <- list()
  
  register_node <- function(x) {
    if (!is.list(x)) {
      return(-as.integer(x)) # Leaf
    }
    
    this_id <- length(nodes)
    # Placeholder
    nodes[[this_id + 1]] <<- list(q = x$q, aye = x$aye, nay = x$nay)
    
    # We must process the current node's children *after* # the children of nodes already in the list to keep indices predictable.
    return(this_id)
  }
  
  # To avoid complexity, we'll use a standard recursive flattening
  # that builds the table by simply tracking the "next available row"
  
  flat_table <- matrix(0, ncol = 3, nrow = 0)
  
  walk <- function(x) {
    if (!is.list(x)) return(-as.integer(x))
    
    this_row_idx <- nrow(flat_table)
    # Reserve the row with a dummy
    flat_table <<- rbind(flat_table, c(as.integer(sub("sp", "", x[[1]])), 0, 0))
    
    # Recurse
    aye_val <- walk(x$aye)
    nay_val <- walk(x$nay)
    
    # Update the reserved row
    # (using this_row_idx + 1 because R is 1-indexed)
    flat_table[this_row_idx + 1, 2] <<- aye_val
    flat_table[this_row_idx + 1, 3] <<- nay_val
    
    return(this_row_idx)
  }
  
  walk(node)
  return(as.vector(t(flat_table)))
}


pecMap <- .EntropyMap(pecScores[pecValid], PAMap(pecSplits))
balMap <- .EntropyMap(balScores[balValid], PAMap(balSplits))


x <- .FlattenMap(pecMap) |> matrix(3) |> t()
head(x)

# Generate the C++ definition
DecisionTreeLine <- function(name, map) {
  flat <- .FlattenMap(map)
  sprintf("inline constexpr int %s7_SCORES[] = {%s};", 
          name, paste(flat, collapse = ","))
}




# tree_vec: The flat integer vector generated by .FlattenMapPure
# present_splits: A vector of integers (e.g., c(13, 15, 48)) 
#                 representing the splits found in the current case.
DebugWalkTree <- function(tree_vec, present_splits) {
  # Convert present_splits to a logical lookup for speed (like our char array)
  lookup <- logical(256)
  lookup[present_splits + 1] <- TRUE # +1 because R is 1-indexed
  
  cursor <- 0
  path <- c() # To track our steps for debugging
  
  while (TRUE) {
    base <- (cursor * 3) + 1 # +1 for R indexing
    
    split_idx <- tree_vec[base]
    aye_target <- tree_vec[base + 1]
    nay_target <- tree_vec[base + 2]
    
    # Log the step
    is_present <- lookup[split_idx + 1]
    path <- c(path, sprintf("Node %d: (sp%d is %s)", cursor, split_idx, is_present))
    
    # Decide next step
    next_node <- if (is_present) aye_target else nay_target
    
    if (next_node < 0) {
      cat("Path taken:\n", paste(path, collapse = " -> "), "\n")
      return(list(result = -next_node, path = path))
    }
    
    cursor <- next_node
  }
}

pecFlat <- .FlattenMap(pecMap)
for (i in seq_len(sum(pecValid))) {
  message(paste0(pecSplits[, i], collapse = "-"))
  expect_equal(DebugWalkTree(pecFlat, pecSplits[, i])$result,
               pecScores[pecValid][i])
}

header_content <- paste0(
  c("// Generated in data-raw/spr-exact-7.R",
    DecisionTreeLine("PEC", pecMap),
    DecisionTreeLine("BAL", balMap))
)

writeLines(header_content, "src/spr/lookup7.h")




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
balDF <- data.frame(key = balPack, score = balScores[balValid])
# pecDF <- pecDF[order(pecDF$key), ]
# balDF <- balDF[order(balDF$key), ]


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
