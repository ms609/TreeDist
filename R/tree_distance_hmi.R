phylo_to_nested <- function(tree) {
  stopifnot(inherits(tree, "phylo"))
  
  # Ensure tree is rooted and binary (ape usually handles this)
  edge <- tree$edge
  nTips <- length(tree$tip.label)
  
  # Build adjacency list
  children <- vector("list", nTips + tree$Nnode)
  for (i in seq_len(nrow(edge))) {
    parent <- edge[i, 1]
    child  <- edge[i, 2]
    children[[parent]] <- c(children[[parent]], child)
  }
  
  # Recursive builder
  build <- function(node) {
    if (node <= nTips) {
      # It's a leaf → return its label
      return(tree$tip.label[node])
    } else {
      # Internal node → return a list of children
      return(lapply(children[[node]], build))
    }
  }
  
  root <- nTips + 1
  build(root)
}
phylo_to_nested_python_like <- function(tree) {
  newick <- ape::write.tree(tree)
  parse_newick <- function(text) {
    text <- gsub("\\s+","", text)
    text <- gsub(";","", text)
    if (!grepl("^\\(", text)) return(text)
    text <- substring(text, 2, nchar(text)-1)
    parts <- character(); cur <- ""; depth <- 0
    for (i in seq_len(nchar(text))) {
      ch <- substr(text,i,i)
      if (ch=="(") {depth<-depth+1; cur<-paste0(cur,ch)}
      else if (ch==")") {depth<-depth-1; cur<-paste0(cur,ch)}
      else if (ch=="," && depth==0) {parts<-c(parts,cur); cur<-""}
      else cur<-paste0(cur,ch)
    }
    if (nchar(cur)>0) parts<-c(parts,cur)
    lapply(parts, parse_newick)
  }
  parse_newick(newick)
}