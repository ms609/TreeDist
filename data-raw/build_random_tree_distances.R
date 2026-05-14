# Build the slim CID/PID reference table shipped with TreeDist.
#
# Source: the full 24-method `randomTreeDistances` array maintained in
# TreeDistData (loaded from the acid-aphid TreeDistData worktree).  We extract
# only the `cid` and `pid` rows, keep `mean` and `sd`, and record the
# replicate count contributing to each cell (per-`n`, taken from the source's
# `n_repls_cid` / `n_repls_pid` attributes -- written by the topup pipeline).
# See `analysis/reference-data-audit.md` in the acid-aphid repo for
# justification.
#
# Pipeline:
#   1. Run `data-raw/enumerate_small_n.R` first (slow, independent) to
#      produce `data-raw/exact_small_n.rds`.
#   2. Run this script (`source('data-raw/build_random_tree_distances.R')`)
#      from the package root: it pulls MC means/sds from the topped-up
#      TreeDistData source, splices in the exact n=4..7 values, and writes
#      `data/randomTreeDistances.rda` via `usethis::use_data()`.

source_path <- normalizePath(file.path(
  "..", "acid-aphid-TreeDistData", "data", "randomTreeDistances.rda"
), mustWork = TRUE)

src_env <- new.env(parent = emptyenv())
load(source_path, envir = src_env)
src <- src_env[["randomTreeDistances"]]

stopifnot(
  is.array(src),
  length(dim(src)) == 3L,
  all(c("cid", "pid") %in% dimnames(src)[[1L]]),
  all(c("mean", "sd") %in% dimnames(src)[[2L]])
)

n_chr <- dimnames(src)[[3L]]
n_int <- as.integer(n_chr)
stopifnot(!anyNA(n_int))

# Per-n replicate counts from the topup pipeline.  Fall back to 1000
# (the original `randomTreeDistances` budget) wherever the attribute is
# missing -- this keeps the script robust to older source builds.
default_reps <- 1000L
attr_or_default <- function(name) {
  attr_val <- attr(src, name, exact = TRUE)
  if (is.null(attr_val)) {
    rep(default_reps, length(n_chr))
  } else {
    as.integer(attr_val[n_chr])
  }
}
n_repls_cid <- attr_or_default("n_repls_cid")
n_repls_pid <- attr_or_default("n_repls_pid")
stopifnot(length(n_repls_cid) == length(n_chr),
          length(n_repls_pid) == length(n_chr),
          !anyNA(n_repls_cid), !anyNA(n_repls_pid))

randomTreeDistances <- data.frame(
  n = n_int,
  cid_mean = unname(src["cid", "mean", ]),
  cid_sd = unname(src["cid", "sd", ]),
  pid_mean = unname(src["pid", "mean", ]),
  pid_sd = unname(src["pid", "sd", ]),
  n_repls_cid = n_repls_cid,
  n_repls_pid = n_repls_pid,
  row.names = NULL,
  stringsAsFactors = FALSE
)

# Splice in exact enumeration results for n = 4..7.  Produced by
# `data-raw/enumerate_small_n.R`; convention is "include-self" mean over all
# N^2 ordered pairs of labelled rooted binary trees (matches the
# independent-draw MC protocol; see header of `enumerate_small_n.R`).
exact_path <- "data-raw/exact_small_n.rds"
if (!file.exists(exact_path)) {
  stop("Missing `", exact_path,
       "`.  Run `data-raw/enumerate_small_n.R` first.")
}
exact <- readRDS(exact_path)
stopifnot(all(c("n", "cid_mean", "cid_sd", "pid_mean", "pid_sd") %in%
              names(exact)),
          all(exact[["n"]] %in% randomTreeDistances[["n"]]))

for (i in seq_len(nrow(exact))) {
  hit <- match(exact[["n"]][[i]], randomTreeDistances[["n"]])
  randomTreeDistances[hit, "cid_mean"] <- exact[["cid_mean"]][[i]]
  randomTreeDistances[hit, "cid_sd"]   <- exact[["cid_sd"]][[i]]
  randomTreeDistances[hit, "pid_mean"] <- exact[["pid_mean"]][[i]]
  randomTreeDistances[hit, "pid_sd"]   <- exact[["pid_sd"]][[i]]
  # Inf flags "exact"; downstream consumers can treat as infinite precision.
  randomTreeDistances[hit, "n_repls_cid"] <- Inf
  randomTreeDistances[hit, "n_repls_pid"] <- Inf
}

stopifnot(
  nrow(randomTreeDistances) == length(n_chr),
  !anyNA(randomTreeDistances[, c("n", "cid_mean", "cid_sd",
                                 "pid_mean", "pid_sd")]),
  all(randomTreeDistances[["cid_mean"]] >= 0),
  all(randomTreeDistances[["cid_mean"]] <= 1),
  all(randomTreeDistances[["pid_mean"]] >= 0),
  all(randomTreeDistances[["pid_mean"]] <= 1)
)

usethis::use_data(randomTreeDistances, internal = FALSE,
                  compress = "xz", overwrite = TRUE)
