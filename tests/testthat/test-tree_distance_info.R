library("TreeTools", quietly = TRUE)

# Independent reference implementations (plain log2), used to avoid tautology.
.Hb <- function(k, n) {
  p <- k / n
  ifelse(p <= 0 | p >= 1, 0, -(p * log2(p) + (1 - p) * log2(1 - p)))
}
.Hjoint <- function(cells, n) {
  p <- cells / n
  p <- p[p > 0]
  -sum(p * log2(p))
}
.VI <- function(a, b, c, d, n) {
  2 * .Hjoint(c(a, b, c, d), n) - .Hb(a + b, n) - .Hb(a + c, n)
}

# Brute force over every ordered composition; every pairing is covered.
.BruteMaxVI <- function(n) {
  best <- -Inf
  for (a in seq_len(n - 3L)) {
    for (b in seq_len(n - a - 2L)) {
      for (cc in seq_len(n - a - b - 1L)) {
        best <- max(best, .VI(a, b, cc, n - a - b - cc, n))
      }
    }
  }
  best
}

# Build a tree pair related by one NNI move, with subtrees of given sizes.
# `Shape` controls each subtree's topology (matters only when normalizing).
.MinShapeNewick <- function(size, labs, split) {
  if (size == 1L) return(labs[[1]])
  i <- split[size]
  paste0("(", .MinShapeNewick(i, labs[seq_len(i)], split), ",",
         .MinShapeNewick(size - i, labs[(i + 1):size], split), ")")
}
.NNIPair <- function(parts, minShape = FALSE) {
  n <- sum(parts)
  split <- if (minShape) .MinSubtreeEntropy(n)[["split"]] else rep(1L, n)
  labs <- paste0("t", seq_len(n))
  ends <- cumsum(parts)
  starts <- c(1L, head(ends, -1L) + 1L)
  frag <- lapply(seq_len(4L), function(j)
    .MinShapeNewick(parts[j], labs[starts[j]:ends[j]], split))
  list(
    a = ape::read.tree(text = paste0("((", frag[[1]], ",", frag[[2]], "),(",
                                      frag[[3]], ",", frag[[4]], "));")),
    b = ape::read.tree(text = paste0("((", frag[[1]], ",", frag[[3]], "),(",
                                      frag[[2]], ",", frag[[4]], "));"))
  )
}

test_that("NNIMaxStep() matches ClusteringInfoDistance for real NNI moves", {
  for (parts in list(c(1, 1, 1, 1), c(2, 1, 1, 1), c(2, 2, 1, 1), c(2, 2, 2, 1),
                     c(3, 3, 2, 2), c(5, 5, 3, 3), c(1, 4, 5, 10), c(4, 3, 2, 9))) {
    pair <- .NNIPair(parts)
    moveDist <- ClusteringInfoDistance(pair[["a"]], pair[["b"]], normalize = FALSE)
    expect_equal(moveDist, .VI(parts[1], parts[2], parts[3], parts[4], sum(parts)),
                 tolerance = 1e-9)
    # The move is a genuine single NNI: exactly one split differs.
    expect_equal(RobinsonFoulds(pair[["a"]], pair[["b"]]), 2)
  }
})

test_that("NNIMaxStep() equals the brute-force maximum", {
  for (n in 4:14) {
    expect_equal(as.numeric(NNIMaxStep(n)), .BruteMaxVI(n),
                 tolerance = 1e-9, info = paste("n =", n))
  }
})

test_that("NNIMaxStep() is exactly two bits iff n is a multiple of four", {
  n <- 4:40
  val <- as.numeric(NNIMaxStep(n))
  expect_equal(val[n %% 4 == 0], rep(2, sum(n %% 4 == 0)))
  expect_true(all(val[n %% 4 != 0] < 2))
  # Ceiling: no move ever exceeds two bits.
  expect_true(all(val <= 2 + 1e-9))
})

test_that("NNIMaxStep() reports the maximizing configuration", {
  m6 <- NNIMaxStep(6)
  expect_equal(attr(m6, "subtrees"), c(1, 1, 2, 2))
  expect_equal(attr(m6, "splits"), c(2, 3))
  # A vector of leaf counts yields a list per attribute, one entry per tree.
  m68 <- NNIMaxStep(c(6, 8))
  expect_length(attr(m68, "subtrees"), 2L)
  expect_equal(attr(m68, "subtrees")[[2]], c(2, 2, 2, 2))
  expect_equal(attr(m68, "splits")[[1]], c(2, 3))
})

test_that("Normalized NNIMaxStep() matches the package", {
  # n >= 16 with optimal subtrees of size >= 4 exercises the entropy-minimizing DP.
  for (n in c(4, 6, 7, 8, 12, 16, 24)) {
    fnVal <- NNIMaxStep(n, normalize = TRUE)
    parts <- attr(fnVal, "subtrees")
    pair <- .NNIPair(parts, minShape = TRUE)
    expect_equal(as.numeric(fnVal),
                 ClusteringInfoDistance(pair[["a"]], pair[["b"]], normalize = TRUE),
                 tolerance = 1e-9, info = paste("n =", n))
  }
  # Four maximally different 4-leaf trees share no clustering information.
  expect_equal(as.numeric(NNIMaxStep(4, normalize = TRUE)), 1)
})

test_that("Entropy-minimizing subtree shape beats a caterpillar in the real package", {
  # For a size-6 subtree the minimum-entropy shape is not the balanced one, so the
  # DP-shaped witness must yield a *larger* normalized distance than a pectinate one.
  parts <- attr(NNIMaxStep(24, normalize = TRUE), "subtrees")
  minPair <- .NNIPair(parts, minShape = TRUE)
  catPair <- .NNIPair(parts, minShape = FALSE)
  expect_gt(ClusteringInfoDistance(minPair[["a"]], minPair[["b"]], TRUE),
            ClusteringInfoDistance(catPair[["a"]], catPair[["b"]], TRUE))
})

test_that("NNIMaxStep() accepts trees and leaf counts alike", {
  expect_equal(as.numeric(NNIMaxStep(BalancedTree(19))),
               as.numeric(NNIMaxStep(19)))
  expect_equal(as.numeric(NNIMaxStep(PectinateTree(19))),
               as.numeric(NNIMaxStep(19)))
  # multiPhylo and list dispatch
  trees <- structure(list(BalancedTree(8), PectinateTree(6)), class = "multiPhylo")
  six <- as.numeric(NNIMaxStep(6))
  expect_equal(as.numeric(NNIMaxStep(trees)), c(2, six))
  lst <- NNIMaxStep(list(BalancedTree(8), 6))
  expect_equal(vapply(lst, as.numeric, double(1)),
               c(2, as.numeric(NNIMaxStep(6))))
})

test_that("NNIMaxStep() retains names", {
  expect_equal(names(NNIMaxStep(c(a = 7))), "a")
  expect_equal(names(NNIMaxStep(c(a = 7, b = 9))), c("a", "b"))
  expect_null(names(NNIMaxStep(7)))
})

test_that("NNIMaxStep() is NA where no NNI move exists", {
  expect_equal(as.numeric(NNIMaxStep(1:4)), c(NA, NA, NA, 2))
  expect_null(attr(NNIMaxStep(3), "subtrees"))
  expect_null(attr(NNIMaxStep(3), "splits"))
})

test_that("NNIMaxStep() agrees with SplitEntropy()", {
  # SplitEntropy()[["Hd"]] is the variation of information between two splits.
  for (parts in list(c(2, 2, 1, 1), c(3, 3, 2, 2), c(5, 5, 3, 3), c(4, 2, 3, 1))) {
    n <- sum(parts)
    a <- parts[1]; b <- parts[2]; cc <- parts[3]; d <- parts[4]
    s1 <- logical(n); s1[seq_len(a + b)] <- TRUE                     # split (a+b | c+d)
    s2 <- logical(n); s2[c(seq_len(a), (a + b) + seq_len(cc))] <- TRUE  # split (a+c | b+d)
    expect_equal(unname(SplitEntropy(s1, s2)[["Hd"]]),
                 .VI(a, b, cc, d, n), tolerance = 1e-9)
  }
})

# Reference R implementations, retained from before the C++ port (#204).
.RefHb <- function(k, n) vapply(k, function(k) Ntropy(k, n - k), double(1))
.RefMinSubtreeEntropy <- function(n) {
  cost <- numeric(n)
  split <- integer(n)
  if (n >= 2L) {
    cost[[2]] <- .RefHb(2L, n)
    split[[2]] <- 1L
  }
  for (s in seq_len(n)[-(1:2)]) {
    i <- seq_len(s %/% 2L)
    daughters <- cost[i] + cost[s - i]
    best <- which.min(daughters)
    cost[s] <- .RefHb(s, n) + daughters[best]
    split[s] <- i[best]
  }
  list(cost = cost, split = split)
}
.RefNNIMaxStep1 <- function(n, normalize) {
  cost <- if (normalize) .RefMinSubtreeEntropy(n)[["cost"]] else NULL
  best <- -Inf
  for (p in seq_len(n %/% 4L)) {
    for (q in p:((n - p) %/% 3L)) {
      r <- q:((n - p - q) %/% 2L)
      s <- n - p - q - r
      h1 <- .RefHb(p + q, n)
      h2 <- .RefHb(p + r, n)
      h3 <- .RefHb(p + s, n)
      keptEntropy <- h1 + h2 + h3 - pmax(h1, h2, h3)
      joint <- vapply(seq_along(r), function(i) Ntropy(p, q, r[[i]], s[[i]]),
                      double(1))
      vi <- 2 * joint - keptEntropy
      value <- if (normalize) {
        vi / (2 * (cost[[p]] + cost[[q]] + cost[r] + cost[s]) + keptEntropy)
      } else {
        vi
      }
      winner <- which.max(value)
      if (value[[winner]] > best) {
        best <- value[[winner]]
        rw <- r[[winner]]
        sw <- s[[winner]]
        bestParts <- c(p, q, rw, sw)
        sizes <- c(p + q, p + rw, p + sw)
        entropies <- c(h1, h2[[winner]], h3[[winner]])
        bestSplits <- sort(sizes[-which.max(entropies)])
      }
    }
  }
  list(value = best, subtrees = bestParts, splits = bestSplits)
}

# Exact ties can be broken differently where the compiler contracts to FMA
# (e.g. ARM64), so compare values to a tolerance and check that the reported
# configuration attains the optimum, rather than demanding identical choices.
test_that(".MinSubtreeEntropy() matches R reference", {
  for (n in c(0:6, 17, 64)) {
    ours <- .MinSubtreeEntropy(n)
    ref <- .RefMinSubtreeEntropy(as.integer(n))
    cost <- ours[["cost"]]
    expect_equal(cost, ref[["cost"]], tolerance = 1e-12)
    split <- ours[["split"]]
    expect_identical(split[seq_len(min(n, 2))], ref[["split"]][seq_len(min(n, 2))])
    for (sz in seq_len(n)[-(1:2)]) {
      i <- seq_len(sz %/% 2L)
      expect_true(split[[sz]] %in% i)
      expect_equal(cost[[split[[sz]]]] + cost[[sz - split[[sz]]]],
                   min(cost[i] + cost[sz - i]), tolerance = 1e-12)
    }
  }
})

test_that("NNIMaxStep() matches R reference", {
  for (n in c(4:21, 40, 63)) {
    cost <- .MinSubtreeEntropy(n)[["cost"]]
    for (norm in c(FALSE, TRUE)) {
      ours <- .NNIMaxStep1(n, norm)
      ref <- .RefNNIMaxStep1(as.integer(n), norm)
      expect_equal(ours[["value"]], ref[["value"]], tolerance = 1e-12)
      
      parts <- ours[["subtrees"]]
      expect_type(parts, "integer")
      expect_equal(sum(parts), n)
      expect_false(is.unsorted(parts))
      sizes <- parts[[1]] + parts[-1]
      h <- .Hb(sizes, n)
      keptH <- sort(h)[1:2]
      splits <- ours[["splits"]]
      expect_type(splits, "integer")
      expect_false(is.unsorted(splits))
      expect_true(all(splits %in% sizes))
      expect_equal(sort(.Hb(splits, n)), keptH, tolerance = 1e-12)
      vi <- 2 * .Hjoint(parts, n) - sum(keptH)
      expect_equal(if (norm) vi / (2 * sum(cost[parts]) + sum(keptH)) else vi,
                   ours[["value"]], tolerance = 1e-12)
    }
  }
})

test_that("Block bounds give exact results at large n", {
  n <- 1100L # Block bounds are used once s / 2 >= 512
  ours <- .MinSubtreeEntropy(n)
  ref <- .RefMinSubtreeEntropy(n)
  cost <- ours[["cost"]]
  expect_equal(cost, ref[["cost"]], tolerance = 1e-12)
  split <- ours[["split"]]
  optimal <- vapply(3:n, function(sz) {
    i <- seq_len(sz %/% 2L)
    abs(cost[[split[[sz]]]] + cost[[sz - split[[sz]]]] -
          min(cost[i] + cost[sz - i])) < 1e-12
  }, logical(1))
  expect_true(all(optimal))
  
  res <- .NNIMaxStep1(n, TRUE)
  # Value from exhaustive search
  expect_equal(res[["value"]], 0.017407858108181639, tolerance = 1e-12)
  parts <- res[["subtrees"]]
  expect_equal(sum(parts), n)
  keptH <- sort(.Hb(parts[[1]] + parts[-1], n))[1:2]
  vi <- 2 * .Hjoint(parts, n) - sum(keptH)
  expect_equal(vi / (2 * sum(cost[parts]) + sum(keptH)), res[["value"]],
               tolerance = 1e-12)
})
