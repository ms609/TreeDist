benchFiles <- list.files("benchmark", "^bench\\-.*\\.R$", full.names = TRUE)
for (benchFile in sample(benchFiles)) {
  system2("Rscript", c("-e", shQuote(paste0("source('", benchFile, "')"))))
}
