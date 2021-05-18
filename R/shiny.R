#' Graphical user interface for mapping distances and analysing
#' tree space
#' 
#' `MapTrees()` launches a 'Shiny' application for the visualization and
#' evaluation of tree spaces.
#' 
#' # Input tab
#' 
#' The input tab allows for the upload of sets of phylogenetic trees from file.
#' Trees at the start or end of a file can be excluded, and the number of trees
#' can be brought down to a manageable number by uniformly subsampling every
#' _n_th tree.  Samples of c. 100 trees can be analysed in seconds;
#' analysis of larger samples will take longer, particularly with slower
#' methods (e.g. quartet distances; Kruskal-1 MDS; large minimum spanning trees).
#' 
#' Different batches can be plotted with different colours / symbols.
#' 
#' If each tree is associated with a property -- for example, the data or method 
#' used to generate it, or its stratigraphic congruence -- a list of properties 
#' for each tree, with one entry per line/row, can be uploaded along with
#' the trees.  Points in tree space can then be styled according to the 
#' corresponding property.
#' 
#' If trees are subsampled (using the 'Sample every' slider), then the values
#' in the tree properties file can also be subsampled accordingly.
#' Unfortunately there is not yet support for multiple point property files;
#' one file will be applied to all trees, in the sequence that they were added
#' to memory.
#' 
#' 
#' # Analysis tab
#' 
#' Select from a suite of distance methods: clustering information and 
#' phylogenetic information are quick and satisfactory; quartet is slow but
#' gives slightly better mappings; path is very fast but may not reflect
#' evolutionary signal very well; and Robinson--Foulds should probably never
#' be used for analysis; it is included for comparison purposes.
#' 
#' Principle components mappings should suffice for most purposes; 
#' Sammon and Kruskal mappings are slower and seldom differ by much,
#' in character or quality, but may emphasize outliers more.
#' 
#' Partitioning around medoids or minimax-linkage hierarchical clustering
#' will typically find a close-to-optimal clustering where one exists;
#' select additional methods for a more exhaustive search.  
#' To avoid redundant calculation, clusterings are only updated when
#' 'recalculate clustering' is clicked, or the 'maximum cluster number' slider
#' is modified; clustering solutions using more than this many clusters are
#' not considered
#' Clusterings with silhouette coefficients < 0.25 are unlikely to represent 
#' genuine structure and are not reported or depicted.
#' 
#' 
#' # Display tab
#' 
#' Up to 15 dimensions can be depicted; the quality of a mapping -- that is, 
#' the faithfulness of mapped distances to true tree-to-tree distances --
#' is quantified by the product of the Trustworthiness and Continuity metrics,
#' which should exceed 0.9 (at least).
#' 
#' An interactive 3D plot can be explored by dragging the mouse and scrolling,
#' but do be careful to check that three dimensions are enough to depict your
#' data accurately.
#' 
#' The minimum spanning tree is the shortest possible line selecting the chosen
#' subsample of trees; if it takes a convoluted zig-zagging route, then the 
#' mapping is doing a poor job of reflecting true tree to tree distances.
#' 
#' Convex hulls are the smallest polygons enclosing all points in each cluster;
#' they are handy for spotting clusters, but their area does not correspond
#' to a genuine quantity, so should not be interpreted.
#' 
#' Tree numbers correspond to the sequence of trees in their original input 
#' file, before subsampling.
#' 
#' Each tree is denoted by a point, whose symbol can be styled according to
#' cluster membership or according to the file that contains the tree,
#' with each click of 'Add to existing' on the input tab constituting a 
#' new batch with a new symbol.
#' 
#' Points can be coloured according to a category -- the cluster or batch to
#' which they belong, or custom data provided in the Point Property File
#' on the input tab -- or continuously, either by the sequence in which they
#' were added to memory, or according to custom data.
#' 
#' 
#' # Exporting tree spaces
#' 
#' A mapping can be saved to PDF or as a PNG bitmap at the size selected.
#' 
#' 
#' # References
#' 
#' A list of references employed when constructing the tree space is populated
#' according to the methods used; it would be appropriate to cite and briefly 
#' discuss these studies in any publication using figures generated using
#' this application.  The application itself can be cited using 
#' Smith (2020, 2021) below.
#'  
#' 
#' 
#' @seealso
#' Full detail of tree space analysis in R is provided in the accompanying
#' [vignette](https://ms609.github.io/TreeDist/articles/treespace.html).
#' 
#' @references 
#' \insertRef{SmithDist}{TreeDist}
#' 
#' \insertRef{SmithSpace}{TreeDist}
#' 
#' @template MRS
#' @family tree space functions
#' 
#' @importFrom shiny runApp
#' @importFrom shinyjs useShinyjs
#' @export
MapTrees <- function() {
  appDir <- system.file("treespace", package = "TreeDist")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing 'TreeDist'.", 
         call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}

#' @rdname MapTrees
#' @export
Project <- function () {
  .Deprecated('MapTrees')
  MapTrees()
}
