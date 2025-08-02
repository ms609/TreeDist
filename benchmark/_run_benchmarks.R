source("benchmark/_init.R")
sapply(list.files("benchmark", "^bench\\-.*\\.R$", full.names = TRUE), source)
