# As spelling::spell_check_test(vignettes = TRUE, error = TRUE), except that
# words found only in NEWS.md, which routinely names new technical terms, are
# reported without failing the check.  Add them to inst/WORDLIST at release.
notCran <- tolower(Sys.getenv("NOT_CRAN")) %in% c("1", "yes", "true")
if (notCran && requireNamespace("spelling", quietly = TRUE)) {
  pkgDir <- list.files("../00_pkg_src", full.names = TRUE)
  if (!length(pkgDir)) {
    checkDir <- dirname(getwd())
    if (grepl("\\.Rcheck$", checkDir)) {
      sourceDir <- sub("\\.Rcheck$", "", checkDir)
      if (file.exists(sourceDir)) pkgDir <- sourceDir
    }
  }
  if (!length(pkgDir) && basename(getwd()) == "tests" &&
      file.exists("../DESCRIPTION")) {
    pkgDir <- dirname(getwd())
  }
  # Package source is not found from covr's temporary install; skip quietly.
  if (length(pkgDir)) {
    results <- spelling::spell_check_package(pkgDir[1], vignettes = TRUE)
    newsOnly <- vapply(results[["found"]],
                       function(where) all(startsWith(where, "NEWS.md:")),
                       logical(1))
    if (any(newsOnly)) {
      cat("Words only in NEWS.md (not an error):",
          paste(results[["word"]][newsOnly], collapse = ", "), "\n")
    }
    if (any(!newsOnly)) {
      stop("Potential spelling errors: ",
           paste(results[["word"]][!newsOnly], collapse = ", "), call. = FALSE)
    }
  }
}
