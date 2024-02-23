library("TreeTools", quietly = TRUE)
devtools::load_all()
set.seed(0)

message(Sys.time(), ": Starting.")

postTrees <- Postorder(as.phylo(0:200, 182))
message(Sys.time(), ": Distances.")
xx <- PathDist(postTrees)

message(Sys.time(), ": phangorn.")
xx <- phangorn::path.dist(postTrees)
message(Sys.time(), ": End.")
