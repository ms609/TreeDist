# TreeDist R Package Development

TreeDist is an R package providing efficient implementations of functions for the comparison of phylogenetic trees. It includes C++ code for performance, comprehensive testing, benchmarking infrastructure, and extensive CI/CD workflows.
Correctness, speed and user-friendliness are priorities.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.


## Context: Object structures

TreeDist handles phylogenetic trees, which are typically presented in one of two formats:

- A `phylo` object is a list containing elements:
  - `"edge": a 2d matrix in which each row identifies the parent vertex (column 1) and child vertex (column 2) of an edge
    in the tree
  - `"Nnode"`: Number of internal nodes
  - `"tip.label"`: Labels for tips with index 1, 2, ... `NTip(tree) == length(tree[["tip.label"]])`
  A phylo object with the attribute "preorder" has edges and internal nodes listed in a strict preorder sequence
  (see `TreeTools::Preorder()`; the attribute "postorder" indicates that edges are numbered in an arbitrary postorder sequence.

- A `Splits` object (from `TreeTools::as.Splits()`) is a raw matrix where each row corresponds to an edge in the tree.
  Each row is named with the integer index of a node associated with the edge (=bipartition split); the bits of the raw
  vector determine which of the two bipartitions each of the `attr(x, "nTip")` leaves (labelled as in the "tip.label" attribute)
  belongs to.  The `TRUE/FALSE` labelling is arbitrary unless `TreeTools::PolarizeSplits` is called.

TreeDist can be considered a "descendant package" of TreeTools, which I also maintain;
TreeTools contains core functionality for tree and split manipulation, designed to be compatible with TreeDist.


## Working Effectively

### Bootstrap and Development Setup
- Install R and development tools:
  ```bash
  sudo apt update && sudo apt install -y r-base r-base-dev build-essential
  sudo apt install -y libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev
  ```
- Install R development packages:
  ```bash
  sudo R -e "install.packages(c('devtools', 'testthat', 'roxygen2', 'lintr'), repos='https://cran.r-project.org/')"
  ```
- Install TreeDist dependencies:
  ```bash
  sudo R -e "install.packages(c('ape', 'bit64', 'lifecycle', 'colorspace', 'fastmatch', 'RCurl', 'R.cache', 'Rdpack', 'stringi', 'PlotTools', 'TreeTools'), repos='https://cran.r-project.org/', dependencies=TRUE)"
  ```

### Building and Checking
- **NEVER CANCEL**: R package builds can take 10-30+ minutes depending on dependencies and system performance. Set timeout to 60+ minutes.
- Build the package:
  ```bash
  R CMD build .
  ```
- **NEVER CANCEL**: Full R CMD check takes 15-45+ minutes. Set timeout to 90+ minutes.
- Check package (full validation):
  ```bash
  R CMD check TreeDist_*.tar.gz
  ```
- Quick check (development, ~2-5 minutes):
  ```bash
  R CMD check --no-build-vignettes --no-manual .
  ```

### Testing
- **NEVER CANCEL**: Test suite has 61+ test files and can take 10-30+ minutes. Set timeout to 60+ minutes.
- Run all tests:
  ```bash
  R -e "devtools::test()"
  ```
- Run tests with testthat directly:
  ```bash
  R -e "testthat::test_dir('tests/testthat')"
  ```
- Run specific test file:
  ```bash
  R -e "testthat::test_file('tests/testthat/test-tree_display.R')"
  ```

### Development Workflow Commands  
- Load package for development:
  ```bash
  R -e "devtools::load_all()"
  ```
- Check code style (~1-2 minutes):
  ```bash
  R -e "lintr::lint_dir('.')"
  ```
- Build documentation:
  ```bash
  R -e "devtools::document()"
  ```
- **NEVER CANCEL**: Build vignettes takes 5-15+ minutes. Set timeout to 30+ minutes.
- Build vignettes:
  ```bash
  R -e "devtools::build_vignettes(install = FALSE)"
  ```

### Key Dependencies
**Critical**: TreeDist requires these packages to build successfully:
- `ape` (>= 5.0) - Phylogenetic analysis package
- `TreeTools` (>= 1.16) - Core tree manipulation (large dependency, 10+ min install)
- `Rdpack` (>= 0.7) - Bibliography and citation support  
- `shinyjs` - Interactive web applications
- `colorspace` - Color space manipulation
- `Rcpp` (>= 1.0.8) - C++ integration (must install before TreeTools)

## Validation

**CRITICAL: Always run the full test suite before proposing changes**
  ```bash
  R -e "devtools::test()"
  ```

- Always run R CMD check for complete validation before finalizing changes.
- ALWAYS run the full test suite when modifying C++ code in src/ directory.
- ALWAYS run lintr to ensure code style compliance before committing.
- For performance-critical changes, run benchmarks in benchmark/ directory:
  ```bash
  R -e "source('benchmark/_run_benchmarks.R')"
  ```
- Memory checking is available but optional (takes significant time):
  ```bash
  R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < memcheck/tests.R
  ```
  **Note**: Memory check scripts available in `memcheck/` directory:
  - `memcheck/tests.R` - Run test suite with valgrind
  - `memcheck/all.R` - Run tests, examples, and build vignettes with memory checking

## Validation Scenarios
After making code changes, validate functionality by testing core phylogenetic tree operations:

### Essential Pre-commit Validation Steps
1. **Install core dependencies first**:
   ```bash
   sudo R -e "install.packages(c('ape', 'colorspace', 'Rdpack', 'shinyjs', 'TreeTools'), repos='https://cran.r-project.org/', dependencies=TRUE)"
   ```

2. **Basic build test** (quick validation):
   ```bash
   R CMD build . --no-build-vignettes
   ```

3. **Run linting** (expects GitHub Actions format):
   ```bash
   R -e "lintr::lint_dir('.')"
   ```

4. **Load package for testing**:
   ```bash
   R -e "devtools::load_all()"
   ```


## Time Expectations & Critical Warnings
- **R startup**: ~0.1 seconds
- **Linting**: 1-3 minutes for full codebase
- **Quick check** (no vignettes/manual): 2-5 minutes  
- **Documentation building**: 2-5 minutes
- **Test suite**: 10-30+ minutes (NEVER CANCEL - set 60+ minute timeout)
- **Full R CMD check**: 15-45+ minutes (NEVER CANCEL - set 90+ minute timeout)
- **Package build**: 10-30+ minutes (NEVER CANCEL - set 60+ minute timeout)
- **Vignette building**: 5-15+ minutes (NEVER CANCEL - set 30+ minute timeout)
- **Benchmarks**: 5-20+ minutes (NEVER CANCEL - set 30+ minute timeout)

## Repository Structure
### Key Directories
- `R/` - R source code (31 R files with 500+ exported functions)
- `src/` - C++ source code requiring compilation (14 C++ files, SystemRequirements: C++17)
- `tests/testthat/` - Test suite (30+ test files)
- `benchmark/` - Performance benchmarking (10+ benchmark files)  
- `man/` - Generated documentation (do not edit manually)
- `vignettes/` - Package tutorials and documentation
- `data/` - Package data files
- `.github/workflows/` - Extensive CI/CD with R-CMD-check, benchmarks, memory checks

### Important Files
- `DESCRIPTION` - Package metadata, dependencies, and system requirements
- `NAMESPACE` - Generated by roxygen2 (do not edit manually)
- `NEWS.md` - Version history (update for user-facing changes)
- `tests/testthat.R` - Test runner entry point
- **Note**: No `.lintr` file exists - linting uses default configuration with GitHub Actions output format

## Common Tasks
### After Making Changes
1. Load and test changes interactively:
   ```bash
   R -e "devtools::load_all(); # test your changes interactively"
   ```
2. Run linting:
   ```bash  
   R -e "lintr::lint_dir('.')"
   ```
3. Run relevant tests:
   ```bash
   R -e "devtools::test()"
   ```
4. For C++ changes, always run full check:
   ```bash
   R CMD check --no-build-vignettes .
   ```

### Code Style Guidelines
- Follow Google's R style guide
- Use camelCase for variable names, TitleCase for function names  
- Use Oxford ending 'ize' (not 'ise') and UK spelling where applicable
- Document functions with roxygen2 comments
- Include test cases for new functionality

### CI Will Fail If
- R CMD check fails
- Tests fail  
- Code style violations (lintr)
- Missing or inadequate documentation
- Missing test coverage for new code

### CI/CD Workflows Available
- **R-CMD-check.yml**: Comprehensive checks on Windows, macOS, Ubuntu across R versions
- **benchmark.yml**: Performance regression testing triggered on PRs
- **memcheck.yml**: Memory checking with valgrind (runs `tests`, `examples`, `vignettes`)
- **ASan.yml**: Address sanitizer checks
- **RcppDeepState.yml**: Deep testing of C++ code
- **pkgdown.yml**: Documentation site generation
- **revdepcheck.yml**: Downstream dependency validation

## Troubleshooting
### Common Build Issues
- **Missing dependencies**: Install system packages first: `sudo apt install -y libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev`
- **Rdpack warnings**: These are normal when `Rdpack` isn't installed but don't prevent building
- **Package won't build**: Install core dependencies: `ape`, `colorspace`, `Rdpack`, `shinyjs`, `TreeTools`
- **TreeTools installation fails**: This is a large dependency - allow 10+ minutes for compilation
- **Tests fail after C++ changes**: rebuild package completely with `R CMD build .`
- **Documentation warnings**: Run `devtools::document()` to regenerate documentation
- **Benchmarks fail**: Performance regressions may need investigation
- **Memory issues**: Use valgrind checking for C++ code validation

### System Requirements Validation
- R version 4.0+ required (tested with R 4.3.3)
- C++17 compiler support required
- Minimum 2GB RAM recommended for building with dependencies
- Allow 60+ minutes for full dependency installation from scratch
