file.pattern <- ".*?(\\w+\\d{4})\\.nex$"
files <- list.files('data-raw', file.pattern, full.names = TRUE)
inapplicable.datasets <- lapply(files, ape::read.nexus.data)
names(inapplicable.datasets) <- inapplicable.names <- gsub(file.pattern, "\\1", files)
inapplicable.phyData <- lapply(inapplicable.datasets, PhyDat)
names(inapplicable.phyData) <- inapplicable.names

devtools::use_data(inapplicable.datasets, inapplicable.phyData, overwrite=TRUE)
