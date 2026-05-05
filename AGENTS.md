# TreeDist — Agent Instructions

## Project Overview

**TreeDist** is an R package (GPL ≥ 3) providing a suite of topological distance metrics
between phylogenetic trees. The mathematical core is implemented in C++17 and exposed to R
via Rcpp. The primary optimization goal is speed: many real analyses compute pairwise
distances over hundreds or thousands of trees, so inner loops must be tight.

Current version: `2.12.0.9000` (development).  
CRAN package page: <https://cran.r-project.org/package=TreeDist>

---

## Repository Layout

```
TreeDist/
├── .AGENTS/              # Protocols for specific activities. Records of work done.
├── src/                  # C++17 source (main optimization target)
├── R/                    # R wrapper layer and pure-R helpers
├── benchmark/            # Microbenchmark scripts (bench package)
├── tests/testthat/       # Unit tests
├── data-raw/             # Scripts that regenerate lookup tables / data
├── vignettes/            # User-facing tutorials
└── inst/                 # Installed extras
```

**Sibling repository** — `../TreeTools` is a companion package that TreeDist links against
at the C++ level (`LinkingTo: TreeTools`). Edits to TreeTools headers (especially
`SplitList.h`) can affect TreeDist performance and can be pushed to CRAN independently
when ready.

---

## C++ Source Files

| File | Size | Purpose |
|------|------|---------|
| `tree_distances.cpp` | 22 KB | Main distance calculations; calls into CostMatrix / LAP |
| `tree_distances.h` | 6 KB | Tree-distance scoring functions (`add_ic_element`, `one_overlap`, `spi_overlap`); includes `lap.h` |
| `lap.h` | 15 KB | **CostMatrix** class; `LapScratch`; `lap()` declarations; cache-aligned storage; `findRowSubmin` hot path |
| `lap.cpp` | 10 KB | Jonker-Volgenant LAPJV linear assignment; extensively hand-optimized; includes only `lap.h` |
| `spr.cpp` | 7 KB | SPR distance approximation |
| `spr_lookup.cpp` | — | SPR lookup-table implementation |
| `nni_distance.cpp` | 16 KB | NNI distance approximations; HybridBuffer allocation |
| `li_diameters.h` | 30 KB | Precomputed NNI diameter lookup tables |
| `information.h` | 6 KB | log₂ / factorial lookup tables (max 8192); cached at startup |
| `binary_entropy_counts.cpp` | — | Binary entropy calculations |
| `day_1985.cpp` | 10 KB | Consensus tree information |
| `hmi.cpp` | 6 KB | Hierarchical Mutual Information |
| `hpart.cpp` | 7 KB | Hierarchical partition structures |
| `reduce_tree.cpp` | 11 KB | Prune trees to common tip set before distance calculation |
| `path_vector.cpp` | 3 KB | Path distance vector |
| `mast.cpp` | 5 KB | Maximum Agreement Subtree |
| `RcppExports.cpp` | 20 KB | Auto-generated Rcpp glue (do not edit by hand) |
| `ints.h` | — | Fixed-width integer typedefs (`splitbit`, `int16`, `int32`, …) |

C++ compilation flags are controlled by `src/Makevars.win` (Windows) / `src/Makevars`.
The package requires **C++17** (`CXX_STD = CXX17`).

---

## Benchmark Infrastructure

When benchmarking, profiling or optimizing, you MUST first read `.AGENTS/protocol/Optimization.md`.

---

## TreeTools Dependency

`../TreeTools` is available locally and editable. It is linked at the C++ level via
`LinkingTo`; the key header consumed by TreeDist is `<TreeTools/SplitList.h>`.

Important constants defined there that affect TreeDist:
- `SL_MAX_TIPS` — maximum number of leaf taxa per tree.
- `SL_MAX_SPLITS` — maximum number of splits.
- `splitbit` — the unsigned integer type used for bitset representation of splits.

If a bottleneck traces back to a TreeTools header (e.g. SplitList layout, bit-width of
`splitbit`), changes can be made in `../TreeTools`, tested locally by re-installing, and
pushed to CRAN when stable.


## Prefer TreeTools over ape

When a function exists in both `ape` and `TreeTools`, **always use the
TreeTools version**. TreeTools wrappers validate input and return
consistently-classed objects.

| ape | TreeTools |
|-----|-----------|
| `reorder(tree, "cladewise")` | `Cladewise(tree)` |
| `reorder(tree, "postorder")` | `Postorder(tree)` |
| `Ntip(tree)` | `NTip(tree)` |
| `keep.tip(tree, tips)` | `KeepTip(tree, tips)` |

Use `ape::` only for functions with no TreeTools equivalent (e.g.
`ape::read.nexus()`, `ape::read.tree()`).


---

## Development Workflow

A task is complete only when R CMD check passes.

Validation on GitHub actions is preferred.
Where there is a strong case to do so, you may test locally:

```r
# Build and reload (from R)
devtools::load_all()          # fast incremental rebuild
devtools::test()              # run testthat suite

# Or from the shell
R CMD build .
R CMD check TreeDist_*.tar.gz

# Run a single benchmark interactively
source("benchmark/_init.R")
source("benchmark/bench-tree-distances.R")
```

## Remote compute dispatch (GHA)

Agents can offload full test suites, R CMD check, and benchmarks to
GitHub Actions instead of running them locally. This frees local CPU
for fast iteration (targeted tests only).

### Helper scripts (in `GitHub/` root)

| Script | Purpose |
|--------|---------|
| `verify-worktree.sh <directory>` | **Mandatory pre-work gate.** Verify worktree branch matches `.worktree-map` |
| `gha-dispatch.sh <workflow> <branch> [input=value ...]` | Trigger a GHA workflow, print run ID |
| `gha-poll.sh <run_id>` | Check run status (exit 0=pass, 1=fail, 2=running) |
| `gha-results.sh <run_id> [output_dir]` | Download artifacts from a completed run |
| `gha-check-pending.sh` | Scan `.gha-pending/` for completed runs ready for pickup |

### Available workflows

| Workflow | Trigger | Purpose |
|----------|---------|---------|
| `agent-check.yml` | `workflow_dispatch` | R CMD check + optional filtered/extended tests (Ubuntu + Windows) |
| `R-CMD-check.yml` | push to main, PRs | Standard CRAN-style check (fires automatically on PRs) |
| `extended-tests.yml` | schedule, `workflow_dispatch` | Tier 3 stress/bench tests |

### Typical usage

```bash
# Push feature branch and dispatch checks
cd TS-<task-name>
git push -u origin feature/<task-name>
cd ..
bash gha-dispatch.sh agent-check.yml feature/<task-name>
# → prints run ID, e.g. 23496860826
```

**Immediately after dispatch, write a pending-file** (see "GHA pending-file
protocol" below). This ensures any agent can pick up the result later:

```bash
RUN_ID=23496860826
cat > .gha-pending/TreeSearch-${RUN_ID}.md << 'EOF'
# GHA Run 23496860826

- **Package:** TreeSearch
- **Branch:** feature/<task-name>
- **Workflow:** agent-check.yml
- **Dispatched by:** Agent <Letter>
- **Dispatched at:** <ISO-8601>
- **Task:** <Letter>-nnn — <description>
- **Worktree:** TS-<task-name>
## Context

<1-3 sentences: what changes are being validated>

## On PASS

<What to do: e.g. "Create PR to cpp-search" or "Report pass on existing PR #N">

## On FAIL

<What to check: e.g. "Timeout logic in ts_driven.cpp line ~443.
Check if elapsed() fires before first replicate completes.">

## Key files changed

- <list of files relevant to diagnosing failures>
EOF
```

### GHA pending-file protocol

**Directory:** `.gha-pending/` (in the `GitHub/` root).

When an agent dispatches a GHA workflow, it writes a context file to
`.gha-pending/<Package>-<run_id>.md`. This file serves two purposes:

1. **Discovery:** Any agent can scan the directory to find completed runs.
2. **Context handoff:** The file contains enough information for a different
   agent to interpret the results and continue the work.

## CPU limits

Max **2 cores** per agent. Use `nThreads = 2L` in tests and benchmarks.
Never use `nThreads = 0L` (auto-detect). Use `-j2` for make.


### Documentation and spelling checks

After any work that adds or modifies roxygen comments, Rd files, NEWS.md, or
vignettes, run:

```r
devtools::check_man()                # regenerates Rd files and checks for issues
spelling::spell_check_package()      # flags potential misspellings
```

Legitimate technical terms, proper nouns, and code identifiers flagged by the
spell checker should be added to `inst/WORDLIST` (one word per line,
alphabetically sorted).  Only fix actual typos in the source.

### Code coverage

Check that existing tests cover all new code.  The GHA test suite uses codecov.
To check locally:

```r
cov <- covr::package_coverage()
covr::report(cov)                        # interactive HTML report
covr::file_coverage(cov, "src/pairwise_distances.cpp")  # per-file summary
```

Aim for full line coverage of new C++ and R code.  If a new code path is not
exercised by the existing test suite, add targeted tests in `tests/testthat/`.

### You are done when:

A completed task: 
- Passes checks
- Has suitable test coverage
- Follows 
