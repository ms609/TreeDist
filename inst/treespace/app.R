suppressPackageStartupMessages({
  library("shiny", exclude = 'runExample')
  library("shinyjs", exclude = 'runExample', warn.conflicts = FALSE)
})
library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
library("TreeDist")

if (!requireNamespace('cluster', quietly = TRUE)) install.packages('cluster')
if (!requireNamespace('protoclust', quietly = TRUE)) install.packages('protoclust')
if (!requireNamespace('MASS', quietly = TRUE)) install.packages('MASS')
if (!requireNamespace('Quartet', quietly = TRUE)) install.packages('Quartet')
if (!requireNamespace('rgl', quietly = TRUE)) install.packages('rgl')
if (!requireNamespace('readxl', quietly = TRUE)) install.packages('readxl')

# Allow large files to be submitted
options(shiny.maxRequestSize = 100 * 1024^2)

palettes <- list("#91aaa7",
                 c("#497a8f", "#c1d0c0"),
                 c("#be83ae", "#2ea7af", "#fbcdcf"),
                 c("#72c5a9", "#b7c5ff", "#dbdabb", "#a28bac"),
                 c("#59c4c0", "#ea9a9a", "#7998a6", "#e9d7a9", "#9c9379"),
                 c("#e8b889", "#67c6eb", "#e5d5c2", "#938fba", "#64b69a", "#779c7b"),
                 c("#c4808f", "#5ba08f", "#f0a693", "#ccd7fe", "#cdb87e", "#c6aae2", "#d2dad8"),
                 c("#d0847f", "#63a5d7", "#d7b981", "#5a9bbb", "#9bb67e", "#dea6d5", "#91967e", "#ca7f96"),
                 c("#8b93a8", "#ccb97e", "#8e9dd7", "#57a384", "#dbb1e7", "#2da7af", "#d68986", "#75d2f9", "#e4d1f0"),
                 c("#dfcf92", "#40b3cb", "#b88a61", "#ecb2e0", "#d6dbbc", "#a28bae", "#edcfeb", "#7498ab", "#b187a0", "#8f939c"),
                 c("#a98f70", "#7be5e2", "#d295c0", "#9ae2bd", "#d3b7f1", "#eca88d", "#8993cd", "#ffc7bb", "#8695a8", "#b3e1df", "#b6878a"),
                 c("#eaa9d3", "#7ac09b", "#fdaca8", "#8ce7e4", "#eed69b", "#70a4d9", "#e8d6ba", "#589bbb", "#959672", "#d0dbd1", "#9b9282", "#d9d9c6"),
                 c("#7498ab", "#e5bd8a", "#7ed8ff", "#de8f8e", "#46bac6", "#ffc0d3", "#57ae96", "#f7cddd", "#50a098", "#b58a6d", "#add49d", "#a18da1", "#cedad9"),
                 c("#8097a4", "#d0dea9", "#a78cc3", "#aee4bf", "#bb82a8", "#5dc9c6", "#b88690", "#26a3b9", "#ad8e6f", "#a4e2ef", "#869a65", "#efcfdd", "#60a089", "#9e927b"),
                 c("#b9aae5", "#bbd69c", "#e2adde", "#77a777", "#f8abc8", "#8ee7ce", "#f2a1a5", "#81bdf1", "#f2bb91", "#b8dcfe", "#aeb276", "#f2cdef", "#e8d6b2", "#8d92b0", "#b7878d"),
                 c("#c3d79b", "#b28cc0", "#64a985", "#e3a7d4", "#2ea2aa", "#e69892", "#85c6f9", "#fbd1a0", "#7696be", "#89996c", "#ddcdff", "#719d89", "#f5cde6", "#b6e0da", "#e8d4cd", "#b5ddfa"),
                 c("#a98d83", "#84e1ff", "#bb8964", "#46b1d1", "#ffbfa5", "#6199c0", "#bbcb8f", "#bf82ab", "#85ddc4", "#eea0ba", "#c1d8ff", "#c3818b", "#c5c6ff", "#999388", "#e8cbff", "#ffb5b6", "#d2dad7"),
                 c("#faccde", "#60a987", "#c6abe4", "#6f9e77", "#c48093", "#a5e5d3", "#cc8f6f", "#499fae", "#d9dca6", "#7796b8", "#bee1ba", "#b4daff", "#919583", "#e2d3e9", "#47a19b", "#ebd4bc", "#7c9993", "#a9e3e0"),
                 c("#739e6e", "#ffbfd9", "#43b6bb", "#e8ad88", "#5e9bce", "#c2af75", "#a8e0fe", "#fad0a8", "#679e8d", "#ffc7b1", "#abe5c0", "#ac8d78", "#c5dddc", "#a48f84", "#cadfb0", "#899694", "#fdcdc1", "#d1dad5", "#dfd8c4"),
                 c("#6e9c93", "#ffb4b3", "#7ec6a2", "#eeccfe", "#cddb9d", "#8a90c5", "#dcb983", "#77bff0", "#f0ab92", "#90ddff", "#f1d3a9", "#b5c2fe", "#c1e1b7", "#7596ba", "#bce1c4", "#a88c96", "#5a9daf", "#b18b80", "#d4d6f3", "#949577"),
                 c("#e7d6bb", "#749ed5", "#f9d29d", "#67b3e2", "#d09773", "#65ccec", "#d38585", "#7fe8ef", "#cf8190", "#94e8cd", "#ae8cc1", "#b3cf95", "#cbc0fc", "#94a66c", "#eeccff", "#ada368", "#e9a6ce", "#48a297", "#ffc1df", "#799c7a", "#facbe0", "#5d9e9a", "#ffc6c1", "#619bb0", "#fccdcb", "#7197bb", "#b1e4c3", "#9390b1", "#c3e0c0", "#a98c90", "#ade3ce", "#9c927d", "#c2dafe", "#869881", "#e6d3dc", "#6e9ba4", "#bde0d0", "#8196a4", "#b2e1df", "#b9deea")
)

badToGood <- rev(c("#1AB958", "#23B956", "#2BB954", "#31B952", "#37B850", "#3CB84E", "#41B84C", "#45B74A", "#49B749", "#4DB747", "#51B645", "#54B643", "#58B641", "#5BB53F", "#5FB53D", "#62B53C", "#65B43A", "#68B438", "#6BB336", "#6DB335", "#70B333", "#73B231", "#76B230", "#78B12E", "#7BB12C", "#7DB02B", "#80B029", "#82AF28", "#85AF26", "#87AE25", "#8AAE23", "#8CAD22", "#8EAD21", "#91AC1F", "#93AC1E", "#95AB1D", "#97AB1C", "#9AAA1B", "#9CAA1A", "#9EA919", "#A0A918", "#A2A818", "#A4A717", "#A6A716", "#A8A616", "#AAA616", "#ACA515", "#AEA415", "#B0A415", "#B2A315", "#B4A315", "#B6A216", "#B8A116", "#B9A117", "#BBA017", "#BD9F18", "#BF9F18", "#C19E19", "#C29D1A", "#C49D1B", "#C69C1C", "#C79B1D", "#C99A1E", "#CB9A1F", "#CC9920", "#CE9822", "#CF9823", "#D19724", "#D29625", "#D49626", "#D59528", "#D79429", "#D8932A", "#D9932C", "#DB922D", "#DC912E", "#DD9130", "#DF9031", "#E08F33", "#E18F34", "#E28E35", "#E38D37", "#E58C38", "#E68C3A", "#E78B3B", "#E88A3D", "#E98A3E", "#EA8940", "#EB8841", "#EC8843", "#ED8744", "#EE8746", "#EE8647", "#EF8549", "#F0854A", "#F1844C", "#F2844D", "#F2834F", "#F38350", "#F48252", "#F48253", "#F58155", "#F58157", "#F68058", "#F6805A", "#F77F5B", "#F77F5D", "#F87E5E"))

Reference <- function (authors, year, title, journal = '',
                       volume = NULL, pages = NULL, doi = NULL,
                       publisher = NULL, editors = NULL) {
  nAuth <- length(authors)
  if (nAuth > 1L) {
    authors <- paste(paste0(authors[-nAuth], collapse = ', '), "&amp;", authors[nAuth])
  }
  nEd <- length(editors)
  if (nEd > 1L) {
    editors <- paste(paste0(editors[-nEd], collapse = ', '), "&amp;", editors[nEd])
  } else if (nEd < 1) {
    editors <- ''
  }
  paste0("<p class='reference'>", authors, " (", year, "). &ldquo;", title,
         "&rdquo;. ",
         if (editors != '') paste0("In: ", editors, ' (eds). ') else '',
         if (journal != "") paste0("<i>", journal, "</i>. ") else "",
         if (is.null(volume)) "" else paste0("<b>", volume, "</b>:"),
         if (is.null(publisher)) "" else paste0(publisher, '. '),
         if (is.null(pages)) "" else paste0(paste0(pages, collapse = "&ndash;"), '. '),
         if (is.null(doi)) "" else paste0(
           "doi:<a href='https://doi.org/", doi, "' title=\"CrossRef\">",
                  doi, "</a>. "), 
         "</p>")
}


Bien2011 <- Reference(
  c("Bien, J.", "Tibshirani, R."),
  title = "Hierarchical clustering with prototypes via minimax linkage",
  year = 2011,
  volume = 106,
  doi = "10.1198/jasa.2011.tm10183",
  pages = c(1075, 1084),
  journal = "Journal of the American Statistical Association")

Day1985 <- Reference(
  title = "Optimal algorithms for comparing trees with labeled leaves",
  author = "Day, W.H.E.", year = 1985,
  volume = 2,
  pages = c(7, 28),
  doi = "10.1007/BF01908061",
  journal = "Journal of Classification")
Estabrook1985 <- Reference(
  c("Estabrook, G.F.", "McMorris, F.R.", "Meacham, C.A."), 1985,
  title = "Comparison of undirected phylogenetic trees based on subtrees of four evolutionary units",
  volume = 34,
  pages = c(193, 200),
  doi = "10.2307/sysbio/34.2.193",
  journal = "Systematic Zoology"
)
Farris1973 <- Reference(title = "On comparing the shapes of taxonomic trees",
                        author = "Farris, J.S.",
                        year = 1973,
                        volume = 22,
                        pages = c(50, 54),
                        doi = "10.2307/2412378",
                        journal = "Systematic Zoology")
Gower1966 <- Reference(  title = "Some distance properties of latent root and vector methods used in multivariate analysis",
                         author = "Gower, J.C.",
                         year = 1966,
                         volume = 53,
                         pages = c(325, 338),
                         doi = "10.2307/2333639",
                         journal = "Biometrika")
Gower1969 <- Reference(
  title = "Minimum spanning trees and single linkage cluster analysis",
  author = c("Gower, J.C.", "Ross, G.J.S."),
  year = 1969, volume = 18, pages = c(54, 64), doi = "10.2307/2346439",
  journal = "Journal of the Royal Statistical Society. Series C (Applied Statistics)")
Kendall2016 <- Reference(
  c('Kendall, M.', 'Colijn, C'), 2016,
  "Mapping phylogenetic trees to reveal distinct patterns of evolution",
  "Molecular Biology and Evolution", 33, c(2735, 2743),
  doi = "10.1093/molbev/msw124")
Kruskal1964 <- Reference(
  title = "Multidimensional scaling by optimizing goodness of fit to a nonmetric hypothesis",
  author = "Kruskal, J.B.", year = 1964, volume = 29, pages = c(1, 27),
  doi = "10.1007/BF02289565", journal = "Psychometrika")
Kaski2003 <- Reference(
  title = "Trustworthiness and metrics in visualizing similarity of gene expression",
  author = c("Kaski, S.", "Nikkil&auml;, J.", "Oja, M.", "Venna, J.",
             "T&ouml;r&ouml;nen, P.", "Castr&eacute;n, E."),
  year = 2003, volume = 4, pages = 48, doi = "10.1186/1471-2105-4-48",
  journal = "BMC Bioinformatics")
Maechler2019 <- Reference(
  title = "cluster: cluster analysis basics and extensions", year = 2019,
  author = c("Maechler, M.", "Rousseeuw, P.", "Struyf, A.", "Hubert, M.", "Hornik, K."),
  journal = "Comprehensive R Archive Network")
Murtagh1983 <- Reference(
  title = "A survey of recent advances in hierarchical clustering algorithms",
  author = "Murtagh, F.", year = 1983, volume = 26, pages = c(354, 359),
  doi = "10.1093/comjnl/26.4.354", journal = "The Computer Journal")
RCoreTeam <- Reference(
  author = "R Core Team", year = 2020,
  title = "R: A language and environment for statistical computing",
  publisher = "R Foundation for Statistical Computing, Vienna, Austria")
Sammon1969 <- Reference(
  title = "A nonlinear mapping for data structure analysis",
  author = "Sammon, J.W.", year = 1969,
  volume = "C-18", pages = c(401, 409),
  doi = "10.1109/T-C.1969.222678", journal = "IEEE Transactions on Computers")
Venables2002 <- Reference(
  title = "Modern Applied Statistics with S. Fourth edition",
  author = c("Venables, W.N.", "Ripley, B.D."),
  publisher = "Springer, New York",
  year = 2002)
Venna2001 <- Reference(
  title = "Neighborhood preservation in nonlinear projection methods: an experimental study",
  author = c("Venna, J.", "Kaski, S."), year = 2001, pages = c(485, 491),
  journal = "Lecture Notes in Computer Science: Artificial Neural Networks&mdash;ICANN 2001",
  editors = c("Dorffner, G.", "Bischof, H.", "Hornik, K."),
  publisher = "Springer, Berlin",
  doi = "10.1007/3-540-44668-0_68")

 
Sand2014 <- Reference(
    author = c("Sand, A.", "Holt, M.K.", "Johansen, J.", "Brodal, G.S.",
               "Mailund, T.", "Pedersen, C.N.S."),
    doi = "10.1093/bioinformatics/btu157",
    journal = "Bioinformatics",
    pages = c(2079, 2080),
    title = "tqDist: a library for computing the quartet and triplet distances between binary or general trees",
    volume = 30,
    year = 2014
)
Smith2020 <- Reference('Smith, M.R.', 2020,
  'Information theoretic Generalized Robinson-Foulds metrics for comparing phylogenetic trees',
  'Bioinformatics', pages = 'In production', doi = "10.1093/bioinformatics/btaa614")
Smith2021 <- Reference('*Smith, M.R.', 2021,
  'The importance of methodology when analyzing landscapes of phylogenetic trees',
  'Submitted MS')
SmithDist <- Reference('Smith, M.R.', 2020,
  'TreeDist: distances between phylogenetic trees',
  doi = '10.5281/zenodo.3528123', 'Comprehensive R Archive Network')
SmithQuartet <- Reference('Smith, M.R.', 2019,
  'Quartet: comparison of phylogenetic trees using quartet and split measures',
  'Comprehensive R Archive Network', doi = "10.5281/zenodo.2536318")
Steel1993 <- Reference(
  c('Steel, M.', 'Penny, D'), 1993,
  "Distributions of tree comparison metrics\u2212some new results",
  "Systematic Biology", 42, c(126, 141),
  doi = "10.1093/sysbio/42.2.126")
Stockham2002 <- Reference(
  author = c('Stockham, C.', 'Wang, L.-S.', 'Warnow, T.'), 2002,
  "Statistically based postprocessing of phylogenetic analysis by clustering",
  "Bioinformatics", 18, c('S285', 'S293'),
  doi = "10.1093/bioinformatics/18.suppl_1.S285")
Ward1963 <- Reference('Ward, J.H.', 1963,
                      'Hierarchical grouping to optimize an objective function',
                      'Journal of the American Statistical Association',
                      58, c(236, 244),
                      doi = "10.1080/01621459.1963.10500845")
Robinson1981 <- Reference(
title = "Comparison of phylogenetic trees",
volume = 53,
doi = "10.1016/0025-5564(81)90043-2",
journal = "Mathematical Biosciences",
author = c("Robinson, D.F.", "Foulds, L.R."),
year = 1981,
pages = c(131, 147)
)


# Define UI for app that draws a histogram ----
ui <- fluidPage(theme = 'treespace.css',
    useShinyjs(),
    column(3, 
      tabsetPanel(
        tabPanel('Input',
          textOutput(outputId = "treesInMemory"),
          fileInput("treefile", "Load trees", placeholder = "No tree file selected"),
          textOutput(outputId = "treesStatus"),
          textOutput(outputId = "keptTrees"),
          sliderInput(inputId = "keepTrees",
                      label = "Trees to analyse:",
                      min = 1,
                      max = 90,
                      value = c(1, 90)),
          sliderInput(inputId = "thinTrees",
                      label = "Sample every:",
                      min = 0,
                      max = 12,
                      round = FALSE,
                      step = 0.01,
                      pre = '2^',
                      value = 0),
          textOutput(outputId = "thinnedTrees"),
          actionButton("replaceTrees", "Replace existing"),
          actionButton("addTrees", "Add batch to existing"),
          
          tags$p("For custom point colours, upload a file with a value for ",
                 "each tree in a separate row: "),
          fileInput('pt.data', 'Point property file',
                    accept = c('.csv', '.txt', '.xls', '.xlsx')),
          checkboxInput('pt.data.subsample', 'Subsample point properties', TRUE),
          tags$p("Subsampling not supported for multiple batches."),
          tagAppendAttributes(textOutput('pt.data.status'), class = 'message'),
          HTML("<p>Do <a href=\"https://github.com/ms609/TreeTools/issues/new?title=Treespace",
               "app:\" title=\"Suggest improvements\">report</a> any bugs or",
               "feature requests to the maintainer (<a",
               "href=\"https://community.dur.ac.uk/martin.smith/\"",
               "title=\"Martin Smith\">Martin R. Smith</a>).</p>"),
        ),
        tabPanel('Analysis',
           selectInput("distance", "Distance method",
                       choices = list("Clustering information" = 'cid',
                                      "Phylogenetic information" = 'pid',
                                      "Quartet (slow)" = 'qd',
                                      "Kendall\u2212Colijn (rooted)" = 'kc',
                                      "Path" = 'path',
                                      "Robinson\u2212Foulds" = 'rf'),
                       selected = 'cid'),
          
          textOutput(outputId = "mappingStatus"),
          fluidRow(plotOutput(outputId = "pcQuality", height = "72px")),
          selectInput("mapping", "Mapping method",
                       choices = list("Princ. comps (classical MDS)" = 'pca',
                                      "Sammon mapping mMDS" = 'nls',
                                      "Kruskal-1 nmMDS (slow)" = 'k'
                                      ),
                       selected = 'pca'),
          
          checkboxGroupInput("clustering", "Clustering methods",
                             choices = list("Partitioning around medoids" = 'pam',
                                            "Hierarchical, minimax linkage" = 'hmm',
                                            "Hierarchical, single linkage" = 'hsi',
                                            "Hierarchical, complete linkage" = 'hco',
                                            "Hierarchical, average linkage" = 'hav',
                                            "Hierarchical, median linkage" = 'hmd',
                                            "Hierarchical, centroid linkage" = 'hct',
                                            "Hierarchical, Ward d\ub2 linkage" = 'hwd',
                                            "K-means" = 'kmn',
                                            "Spectral" = "spec"),
                             selected = c('pam', 'hmm')),
          sliderInput(inputId = "maxClusters",
                      label = "Maximum cluster number:",
                      min = 2,
                      max = 20,
                      value = 8),
          actionButton("calcClusters", "Recalculate clustering"),
          fluidRow(plotOutput(outputId = "clustQuality", height = "80px")),
        ),
        tabPanel('Display',
          sliderInput(inputId = "dims",
                      label = "Dimensions to plot:",
                      min = 2, max = 15, value = 5),
          plotOutput('howManyDims', width = '100%', height = '180px'),
          
          sliderInput(inputId = "mst",
                      label = "Minimum spanning tree size",
                      min = 0,
                      max = 90,
                      value = 90),
          checkboxGroupInput("display", "Display settings",
                             list(
                               "Clusters' convex hulls" = 'hulls',
                               'Cluster consensus trees' = 'cons',
                               'Tree numbers' = 'labelTrees',
                               'Interactive 3D plot' = 'show3d'
                             ),
                             character(0)),
          sliderInput('pt.cex', 'Point size',
                      min = 0, max = 4, value = 1, 
                      step = 0.01),
          sliderInput('pt.opacity', 'Point opacity',
                      min = 0, max = 255, value = 196,
                      step = 1),
          selectInput('pt.pch', 'Point symbols:',
                      list(
                        'Solid circle' = '16',
                        'Open ring' = '1',
                        'One per cluster' = 'clust',
                        'One per batch' = 'entry'
                      ), '16'),
          selectInput('pt.col', 'Colour points by:',
                      list('Cluster membership' = 'clust', 
                           'Batch' = 'entry',
                           'Sequence within batches' = 'seq',
                           'Distance from focal tree' = 'dist',
                           'Custom categorical' = 'discrete',
                           'Custom continuous' = 'continuous'), 'clust'),
          tagAppendAttributes(textOutput('pt.col.status'), class = 'message'),
          plotOutput('pt.col.scale', height = "36px"),
          hidden(sliderInput('which.tree', "Focal tree index",
                             min = 1, max = 90, value = 1, step = 1)),
        )
      ),
    ),
    
    column(9,
      fluidRow(id = 'plotConfig',
        tags$span("Plot size:", id = 'plotSizeSpan'),
        sliderInput(inputId = "plotSize", label = NULL,
                    width = '200px',
                    min = 100, max = 2000,
                    post = 'px', value = 600),
        tags$span("Save as: "),
        downloadButton('savePdf', 'PDF'),
        downloadButton('savePng', 'PNG'),
        downloadButton('saveCons', 'Cluster consensus trees'),
      ),
      fluidRow(
        plotOutput(outputId = "distPlot", height = "0px"),
        rgl::rglwidgetOutput(outputId = "threeDPlot",
                             width = "600px", height = "600px"),
        hidden(plotOutput('clustCons', height = "200px")),
        htmlOutput('references'),
      ),
  )
)

server <- function(input, output, session) {
  #treeNumbers <- c(1:220)
  treeNumbers <- c(1:50, 401:440)
  
  r <- reactiveValues(allTrees = as.phylo(treeNumbers, 8),
                      entryOrder = rep(1, length(treeNumbers)))
  
  ##############################################################################
  # Reset cache
  ##############################################################################
  ClearClusters <- function () {
    r$clust_cid = NULL
    r$clust_pid = NULL
    r$clust_qd = NULL
    r$clust_kc = NULL
    r$clust_path = NULL
    r$clust_rf = NULL
  }
  
  ClearMST <- function () {
    r$mst_cid = NULL
    r$mst_pid = NULL
    r$mst_qd = NULL
    r$mst_kc = NULL
    r$mst_path = NULL
    r$mst_rf = NULL
  }
  
  ClearDistances <- function () {
    r$dist_cid = NULL
    r$dist_pid = NULL
    r$dist_qd = NULL
    r$dist_kc = NULL
    r$dist_path = NULL
    r$dist_rf = NULL
    
    r$proj_pca_cid = NULL
    r$proj_pca_pid = NULL
    r$proj_pca_qd = NULL
    r$proj_pca_kc = NULL
    r$proj_pca_path = NULL
    r$proj_pca_rf = NULL
    
    r$proj_k_cid = NULL
    r$proj_k_pid = NULL
    r$proj_k_qd = NULL
    r$proj_k_kc = NULL
    r$proj_k_path = NULL
    r$proj_k_rf = NULL
    
    r$proj_nls_cid = NULL
    r$proj_nls_pid = NULL
    r$proj_nls_qd = NULL
    r$proj_nls_kc = NULL
    r$proj_nls_path = NULL
    r$proj_nls_rf = NULL
    
    ClearClusters()
    
    ClearMST()
  }
  
  ##############################################################################
  # Load trees
  ##############################################################################
  treeFile <- reactive({
    fileInput <- input$treefile
    if (is.null(input)) return("No tree file selected.")
    tmpFile <- fileInput$datapath
    if (is.null(tmpFile)) return ("No trees found.")
    if (length(grep('#NEXUS', toupper(readLines(tmpFile)[1]),
                    fixed = TRUE)) > 0) {
      ret <- read.nexus(tmpFile)
    } else {
      ret <- ReadTntTree(tmpFile)
      if (length(ret) == 0) ret <- read.tree(tmpFile)
    }
    
    if (!inherits(ret, c('phylo', 'multiPhylo'))) {
      return("Could not read trees from file")
    }
    nTrees <- length(ret)
    if (nTrees) {
      updateSliderInput(session, 'keepTrees',
                        value = c(0, nTrees),
                        max = nTrees)
      updateSliderInput(session, 'thinTrees',
                        value = log2(nTrees / 100),
                        max = max(0, signif(log2(nTrees) - 1L, 3)))
    }
    ret
  })
  
  output$treesStatus <- renderText({
    if (is.character(treeFile())) return(treeFile())
    paste(length(treeFile()), "trees in file.")
  })
  
  treesLoaded <- reactive({
    !is.character(treeFile())
  })
  
  keptRange <- reactive({
    if (treesLoaded()) {
      nTrees <- length(treeFile())
      input$keepTrees
    } else {
      c(1, length(r$allTrees))
    }
  })
  
  output$keptTrees <- renderText({
    nTrees <- length(treeFile())
    if (treesLoaded()) {
      paste0("Keeping trees ", keptRange()[1], " to ", keptRange()[2], ".")
    } else {
      "Using example trees."
    }
  })
  
  output$thinnedTrees <- renderText({
    if (input$thinTrees == 1) "" else {
      nTrees <- length(thinnedTrees())
      if (nTrees > 1500) {
        paste("Selecting", nTrees, "trees.  I hope you're feeling patient.")
      } else if (nTrees > 500) {
        paste("Selecting", nTrees, "trees.  May be slow.")
      } else {
        paste("Selecting", nTrees, "trees.")
      }
    }
  })
  
  thinnedTrees <- reactive({
    as.integer(seq(keptRange()[1], keptRange()[2], by = 2^input$thinTrees))
  })
  
  TreeNumberUpdate <- function () {
    ClearDistances()
    nTrees <- length(r$allTrees)
    updateSliderInput(session, 'mst', value = min(500, nTrees), max = nTrees)
    updateSliderInput(session, 'which.tree', max = nTrees)
    updateSliderInput(session, 'dims',
                      max = nProjDim(),
                      value = min(input$dims, nProjDim()))
  }
  
  observeEvent(input$addTrees, {
    newTrees <- treeFile()[thinnedTrees()]
    if (length(r$allTrees) == 0) {
      r$entryOrder <- rep(1L, length(newTrees))
      r$allTrees <- c(newTrees)
    } else {
      r$entryOrder <- c(r$entryOrder,
                        rep(1L + max(r$entryOrder), length(newTrees)))
      r$allTrees <- c(r$allTrees, newTrees)
    }
    TreeNumberUpdate()
  })
  
  observeEvent(input$replaceTrees, {
    r$allTrees <- c(treeFile()[thinnedTrees()])
    r$entryOrder <- rep(1L, length(r$allTrees))
    TreeNumberUpdate()
  })
  
  output$treesInMemory <- renderText({
    paste(length(r$allTrees), "trees in memory.")
  })
  
  nLeaves <- reactive(max(NTip(r$allTrees)))
  
  ##############################################################################
  # Calculate distances
  ##############################################################################

  chosenDistance <- function() {
    switch(input$distance, 
           cid = "Clustering information distance",
           pid = "Phylogenetic information distance",
           qd = "Quartet divergence",
           kc = "Kendall\u2212Colijn distance",
           path = "Path distance",
           rf = "Robinson\u2212Foulds distance")
  }
  
  distances <- reactive({
    if (length(r$allTrees) > 1L) {
      withProgress(
        message = 'Calculating distances', value = 0.99,
        switch(input$distance,
               'cid' = {
                 if (is.null(r$dist_cid)) r$dist_cid = ClusteringInfoDistance(r$allTrees)
                 r$dist_cid
               },
               'pid' = {
                 if (is.null(r$dist_pid)) r$dist_pid = PhylogeneticInfoDistance(r$allTrees)
                 r$dist_pid
               },
               'qd' = {
                 if (is.null(r$dist_qd)) r$dist_qd = as.dist(Quartet::QuartetDivergence(
                   Quartet::ManyToManyQuartetAgreement(r$allTrees), similarity = FALSE))
                 r$dist_qd
               },
               'kc' = {
                 if (is.null(r$dist_kc)) r$dist_kc = KendallColijn(r$allTrees)
                 r$dist_kc
               },
               'path' = {
                 if (is.null(r$dist_path)) r$dist_path = PathDist(r$allTrees)
                 r$dist_path
               },
               'rf' = {
                 if (is.null(r$dist_rf)) r$dist_rf = RobinsonFoulds(r$allTrees)
                 r$dist_rf
               }
         )
      )
    } else {
      matrix(0, 0, 0)
    }
    
  })
  
  ##############################################################################
  # Mapping
  ##############################################################################
  maxProjDim <- reactive({
    min(15L, length(r$allTrees) - 1L)
  })
  
  nProjDim <- reactive({
    dim(mapping())[2]
  })
  
  dims <- debounce(reactive({
    if (mode3D()) 3L else {
      min(input$dims, maxProjDim())
    }
  }), 400)
  
  mapping <- reactive({
    if (maxProjDim() > 1L) {
      proj_id <- paste0('proj_', input$mapping, '_', input$distance)
      if (is.null(r[[proj_id]])) {
        withProgress(
          message = 'Mapping distances',
          value = 0.99,
          proj <- switch(input$mapping,
                         'pca' = cmdscale(distances(), k = maxProjDim()),
                         'k' = MASS::isoMDS(distances(), k = maxProjDim())$points,
                         'nls' = MASS::sammon(distances(), k = maxProjDim())$points
                         )
        )
        r[[proj_id]] <- proj
      }
      
      # Return:
      r[[proj_id]]
    } else {
      matrix(0, 0, 0)
    }
  })
  
  projQual <- reactive({
    withProgress(message = "Estimating mapping quality", {
      vapply(seq_len(nProjDim()), function (k) {
        incProgress(1 / nProjDim())
        MappingQuality(distances(), dist(mapping()[, seq_len(k)]), 10)
      }, numeric(4))
    })
  })
  
  LogScore <- function (x) {
    (-(log10(1 - x + 1e-2))) / 2
  }
  
  output$pcQuality <- renderPlot({
    par(mar = c(2, 0, 0, 0))
    nStop <- length(badToGood)
    
    plot(seq(0, 1, length.out = nStop), rep(0, nStop),
         pch = 15, col = badToGood,
         ylim = c(-1.5, 2.5),
         ann = FALSE, axes = FALSE)
    
    logScore <- LogScore(projQual()['TxC', dims()])
    lines(rep(logScore, 2), c(-1, 1), lty = 3)
    
    tickPos <- c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0)
    ticks <- LogScore(tickPos)
    
    axis(1, at = ticks, labels = NA, line = 0)
    axis(1, tick = FALSE, at = ticks, labels = tickPos, line = 0)
    axis(1, line = -1, tick = FALSE,
         at = ticks[-1] - ((ticks[-1] - ticks[-length(ticks)]) / 2),
         labels = c('', 'dire', '', "ok", "gd", "excellent"))
    axis(3, at = 0.5, tick = FALSE, line = -2, 
         paste0(dims(), 'D mapping quality (trustw. \ud7 contin.):'))
  })
  
  
  output$howManyDims <- renderPlot({
    par(mar = c(2.5, 2.5, 0, 0), xpd = NA, mgp = c(1.5, 0.5, 0))
    txc <- projQual()['TxC', ]
    nStop <- length(badToGood)
    
    plot(txc, type = 'n', ylim = c(min(txc, 0.5), 1),
         frame.plot = FALSE, axes = FALSE,
         xlab = 'Dimension', ylab = 'Trustw. \uD7 Contin.')
    par(xpd = FALSE)
    axis(1, 1:14)
    axis(2)
    tickPos <- c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0)
    mids <- c(0.6, 0.75, 0.85, 0.925)
    text(rep(15, 4), mids, pos = 2, cex = 0.8,
         col = badToGood[nStop * LogScore(mids)],
         c('Essentially random', 'Dangerous', "Usable", "Good"))
    text(1, 0.975, pos = 4, "Excellent", cex = 0.8, 
         col = badToGood[LogScore(0.975) * nStop])
    for (i in tickPos[-1]) {
      abline(h = i, lty = 3, col = badToGood[LogScore(i) * nStop])
    }
    points(txc, type = 'b')
    txcNow <- txc[dims()]
    
    points(dims(), txcNow, pch = 16, col = badToGood[LogScore(txcNow) * nStop],
           cex = 1.6)
  })
  
  ##############################################################################
  # Clusterings
  ##############################################################################
  maxClust <- reactive(min(input$maxClusters, length(r$allTrees) - 1L))
  clusterings <- reactive({
    maxCluster <- input$maxClusters
    if (!is.null(r$clust_max) && maxCluster != r$clust_max) ClearClusters()
    clust_id <- paste0('clust_', input$distance)
    
    if (is.null(r[[clust_id]])) {
      if (maxClust() > 1) {
        possibleClusters <- 2:maxClust()
        
        hSil <- pamSil <- havSil <- kSil <- specSil <-
          hsiSil <- hcoSil <- hmdSil <- hctSil <- hwdSil <- -99
        dists <- distances()
        
        nMethodsChecked <- length(input$clustering)
        methInc <- 1 / nMethodsChecked
        nK <- length(possibleClusters)
        kInc <- 1 / (nMethodsChecked * nK)
        
        withProgress(message = "Clustering", {
          if ('pam' %in% input$clustering) {
            pamClusters <- lapply(possibleClusters, function (k) {
              incProgress(kInc, detail = 'PAM clustering')
              cluster::pam(dists, k = k)
            })
            pamSils <- vapply(pamClusters, function (pamCluster) {
              incProgress(kInc, detail = 'PAM silhouettes')
              mean(cluster::silhouette(pamCluster)[, 3])
            }, double(1))
            
            bestPam <- which.max(pamSils)
            pamSil <- pamSils[bestPam]
            pamCluster <- pamClusters[[bestPam]]$cluster
          }
          
          if ('hmm' %in% input$clustering) {
            incProgress(methInc / 2, detail = 'minimax clustering')
            hTree <- protoclust::protoclust(dists)
            hClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
            hSils <- vapply(hClusters, function (hCluster) {
              incProgress(kInc / 2, detail = 'minimax silhouettes')
              mean(cluster::silhouette(hCluster, dists)[, 3])
            }, double(1))
            bestH <- which.max(hSils)
            hSil <- hSils[bestH]
            hCluster <- hClusters[[bestH]]
          }
          
          if ('hwd' %in% input$clustering) {
            incProgress(methInc / 2, detail = 'Ward D\ub2 clustering')
            hTree <- stats::hclust(dists, method = 'ward.D2')
            hwdClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
            hwdSils <- vapply(hwdClusters, function (hCluster) {
              incProgress(kInc / 2, detail = 'Ward D\ub2 silhouettes')
              mean(cluster::silhouette(hCluster, dists)[, 3])
            }, double(1))
            bestHwd <- which.max(hwdSils)
            hwdSil <- hwdSils[bestHwd]
            hwdCluster <- hwdClusters[[bestHwd]]
          }
          
          if ('hsi' %in% input$clustering) {
            incProgress(methInc / 2, detail = 'single clustering')
            hTree <- stats::hclust(dists, method = 'single')
            hsiClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
            hsiSils <- vapply(hsiClusters, function (hCluster) {
              incProgress(kInc / 2, detail = 'single silhouettes')
              mean(cluster::silhouette(hCluster, dists)[, 3])
            }, double(1))
            bestHsi <- which.max(hsiSils)
            hsiSil <- hsiSils[bestHsi]
            hsiCluster <- hsiClusters[[bestHsi]]
          }
          
          if ('hco' %in% input$clustering) {
            incProgress(methInc / 2, detail = 'complete clustering')
            hTree <- stats::hclust(dists, method = 'complete')
            hcoClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
            hcoSils <- vapply(hcoClusters, function (hCluster) {
              incProgress(kInc / 2, detail = 'complete silhouettes')
              mean(cluster::silhouette(hCluster, dists)[, 3])
            }, double(1))
            bestHco <- which.max(hcoSils)
            hcoSil <- hcoSils[bestHco]
            hcoCluster <- hcoClusters[[bestHco]]
          }
          
          
          if ('hav' %in% input$clustering) {
            incProgress(methInc / 2, detail = 'average clustering')
            hTree <- stats::hclust(dists, method = 'average')
            havClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
            havSils <- vapply(havClusters, function (hCluster) {
              incProgress(kInc / 2, detail = 'average silhouettes')
              mean(cluster::silhouette(hCluster, dists)[, 3])
            }, double(1))
            bestHav <- which.max(havSils)
            havSil <- havSils[bestHav]
            havCluster <- havClusters[[bestHav]]
          }
          
          
          if ('hmd' %in% input$clustering) {
            incProgress(methInc / 2, detail = 'median clustering')
            hTree <- stats::hclust(dists, method = 'median')
            hmdClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
            hmdSils <- vapply(hmdClusters, function (hCluster) {
              incProgress(kInc / 2, detail = 'median silhouettes')
              mean(cluster::silhouette(hCluster, dists)[, 3])
            }, double(1))
            bestHmd <- which.max(hmdSils)
            hmdSil <- hmdSils[bestHmd]
            hmdCluster <- hmdClusters[[bestHmd]]
          }
          
          
          if ('hct' %in% input$clustering) {
            incProgress(methInc / 2, detail = 'centroid clustering')
            hTree <- stats::hclust(dists ^ 2, method = 'centroid')
            hctClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
            hctSils <- vapply(hctClusters, function (hCluster) {
              incProgress(kInc / 2, detail = 'centroid silhouettes')
              mean(cluster::silhouette(hCluster, dists)[, 3])
            }, double(1))
            bestHct <- which.max(hctSils)
            hctSil <- hctSils[bestHct]
            hctCluster <- hctClusters[[bestHct]]
          }
          
          if ('kmn' %in% input$clustering) {
            incProgress(methInc / 2, detail = 'K-means clustering')
            kClusters <- lapply(possibleClusters, function (k) kmeans(dists, k))
            kSils <- vapply(kClusters, function (kCluster) {
              incProgress(kInc / 2, detail = 'K-means silhouettes')
              mean(cluster::silhouette(kCluster$cluster, dists)[, 3])
            }, double(1))
            
            bestK <- which.max(kSils)
            kSil <- kSils[bestK]
            kCluster <- kClusters[[bestK]]$cluster
          }
          
          if ('spec' %in% input$clustering) {
            spectralEigens <- SpectralEigens(dists,
                                             nn = min(ncol(as.matrix(dists)) - 1L, 10),
                                             nEig = 3L)
            specClusters <- lapply(possibleClusters, function (k) {
              incProgress(kInc / 2, detail = 'spectral clustering')
              cluster::pam(spectralEigens, k = k)
            })
            specSils <- vapply(specClusters, function (cluster) {
              incProgress(kInc / 2, detail = 'spectral silhouettes')
              mean(cluster::silhouette(cluster$cluster, dists)[, 3])
            }, double(1))
            
            bestSpec <- which.max(specSils)
            specSil <- specSils[bestSpec]
            specCluster <- specClusters[[bestSpec]]$cluster
          }
          bestCluster <- c('none', 'pam', 'hmm', 'hsi', 'hco', 'hav', 'hmd', 'hct',
                           'hwd', 'kmn', 'spec')[
                             which.max(c(0, pamSil, hSil, hsiSil, hcoSil, havSil, hmdSil, hctSil,
                                         hwdSil, kSil, specSil))]
          
        })
      } else {
        bestCluster <- 'none'
      }
      
      r[[clust_id]] <- list(method = switch(bestCluster,
                                            pam = 'part. around medoids',
                                            hmm = 'minimax linkage',
                                            hsi = "single linkage",
                                            hco = "complete linkage",
                                            hav = "average linkage",
                                            hmd = "median linkage",
                                            hct = "centroid linkage",
                                            hwd = 'Ward d\ub2 linkage',
                                            kmn = 'K-means',
                                            spec = 'spectral',
                                            'no attempt to find any'),
                            n = 1 + switch(bestCluster,
                                           pam = bestPam,
                                           hmm = bestH,
                                           hsi = bestHsi,
                                           hco = bestHco,
                                           hav = bestHav,
                                           hmd = bestHmd,
                                           hct = bestHct,
                                           hwd = bestHwd,
                                           har = bestHav, kmn = bestK,
                                           spec = bestSpec,
                                           0),
                            sil = switch(bestCluster,
                                         pam = pamSil,
                                         hmm = hSil,
                                         hsi = hsiSil,
                                         hco = hcoSil,
                                         hav = havSil,
                                         hmd = hmdSil,
                                         hct = hctSil,
                                         hwd = hwdSil,
                                         kmn = kSil, spec = specSil, 
                                         0), 
                            cluster = switch(bestCluster,
                                             pam = pamCluster,
                                             hmm = hCluster,
                                             hsi = hsiCluster,
                                             hco = hcoCluster,
                                             hav = havCluster,
                                             hmd = hmdCluster,
                                             hct = hctCluster,
                                             hwd = hwdCluster,
                                             kmn = kCluster,
                                             spec = specCluster,
                                             1)
      )
      r$clust_max <- maxCluster
    }
    r[[clust_id]]
  })
  
  observeEvent(input$calcClusters, {ClearClusters()})
  
  mstSize <- debounce(reactive(input$mst), 100)
  
  mstEnds <- reactive({
    
    if (!is.null(r$mst_size) && mstSize() != r$mst_size) ClearMST()
    mst_id <- paste0('mst_', input$distance)
    if (is.null(r[[mst_id]])) {
      dist <- as.matrix(distances())
      nRows <- dim(dist)[1]
      withProgress(message = 'Calculating MST', {
        selection <- unique(round(seq.int(1, nRows, length.out = max(2L, mstSize()))))
        edges <- MSTEdges(dist[selection, selection])
        edges[, 1] <- selection[edges[, 1]]
        edges[, 2] <- selection[edges[, 2]]
        r[[mst_id]] <- edges
      })
    }
    
    r$mst_size <- mstSize()
    r[[mst_id]]
  })
  
  output$clustQuality <- renderPlot({
    par(mar = c(2, 0.5, 0, 0.5), xpd = NA, mgp = c(2, 1, 0))
    cl <- clusterings()
    clust_id <- paste0('clust_', input$distance)
    sil <- r[[clust_id]]$sil
    if (length(sil) == 0) sil <- -0.5
    nStop <- 400
    range <- c(0.5, 1)
    negScale <- hcl.colors(nStop, 'plasma')[seq(range[1] * nStop, 1,
                                            length.out = nStop * range[1])]
    posScale <- hcl.colors(nStop, 'viridis')
    
    plot(seq(-range[1], range[2], length.out = nStop * sum(range)),
         rep(0, nStop * sum(range)),
         pch = 15, col = c(negScale, posScale),
         ylim = c(-2, 2),
         ann = FALSE, axes = FALSE)
    lines(c(sil, sil), c(-1, 1), lty = 3)
    
    ticks <- c(-0.5, 0, 0.25, 0.5, 0.7, 1)
    axis(1, at = ticks, line = -1)
    axis(1, tick = FALSE, at = ticks[-1] - ((ticks[-1] - ticks[-6]) / 2),
         labels = c("Structure:", "none", "weak", "ok", "strong"),
         line = -2)
    
    axis(1, tick = FALSE, line = 0, at = 0.25,
         labels = paste0("Silhouette coefficient (", round(sil, 3), ")"))
    
    
    axis(3, tick = FALSE, line = -2, at = 0.25,
         labels = if (sil > 0.25) 
           paste0(cl$n, " clusters found with ", cl$method) else 
             paste0("No meaningful clusters found"))
    
  })
  
  observeEvent(input$display, {
    if ('cons' %in% input$display) {
      show('clustCons')
    } else {
      hide('clustCons')
    }
  }, ignoreInit = TRUE, ignoreNULL = FALSE)
  
  consRows <- reactive({
    cl <- clusterings()
    if (cl$sil > 0.25) ceiling(cl$n / 3) else 1L
  })
  
  consSize <- reactive({
    nLeaves() * 12 * consRows()
  })
  
  PlotClusterCons <- function () {
    cl <- clusterings()
    par(mar = c(0.2, 0, 0.2, 0), xpd = NA)
    par(cex = 0.9)
    if (cl$sil > 0.25) {
      par(mfrow = c(consRows(), ceiling(cl$n / consRows())))
      for (i in seq_len(cl$n)) {
        col <- palettes[[min(length(palettes), cl$n)]][i]
        tr <- ape::consensus(r$allTrees[cl$cluster == i])
        tr$edge.length <- rep(1, dim(tr$edge)[1])
        plot(tr, edge.width = 2, font = 1, cex = 1,
             edge.color = col, tip.color = col)
      }
    } else {
      tr <- ape::consensus(r$allTrees)
      tr$edge.length <- rep(1, dim(tr$edge)[1])
      plot(tr,edge.width = 2, font = 1, cex = 1,
           edge.color = palettes[[1]],
           tip.color = palettes[[1]])
    }
  }
  
  output$clustCons <- renderPlot({
    if ('cons' %in% input$display) {
      PlotClusterCons()
    }
  }, height = consSize)
  
  ##############################################################################
  # Plot settings: point style
  ##############################################################################
  PointDataStatus <- function (...) {
    msg <- paste0(...)
    output$pt.data.status <- renderText(msg)
    output$pt.col.status <- renderText(msg)
  }
  
  pointData <- reactive({
    fp <- input$pt.data$datapath
    PointDataStatus("")
    extension <- if(is.null(fp)) '' else substr(fp, nchar(fp) - 3, nchar(fp))
    ret <- switch(extension,
                  '.csv' = read.csv(fp),
                  '.txt' = read.table(fp),
                  '.xls' = readxl::read_excel(fp),
                  'xlsx' = readxl::read_excel(fp),
                  {
                    PointDataStatus("Unrecognized file extension.")
                    matrix(0)
                  }
    )
    
    if (input$pt.data.subsample) ret[thinnedTrees(), 1] else ret[, 1]
  })
  
  ContinuousPtCol <- function (dat, bigDark = FALSE) {
    show('pt.col.scale')
    scale <- substr(hcl.colors(256, 'plasma'), 1, 7)
    if (bigDark) scale <- rev(scale)
    output$pt.col.scale <- renderPlot({
      par(mar = c(1, 1, 0, 1))
      plot(1:256, rep(0, 256), pch = 15, col = scale,
           ylim = c(-1, 0), ann = FALSE, axes = FALSE)
      axis(1, at = c(1, 256), labels = signif(range(dat), 4), line = -1)
    })
    scale[cut(dat, 256)]
  }
  
  observeEvent(input$pt.col, {
    if (input$pt.col == 'dist') show('which.tree') else hide('which.tree')
  })
  
  pointCols <- reactive({
    switch(input$pt.col,
           'clust' = {
             hide('pt.col.scale')
             cl <- clusterings()
             if (cl$sil > 0.25) {
               palettes[[min(length(palettes), cl$n)]][cl$cluster]
             } else palettes[[1]]
           },
           'entry' = {
             hide('pt.col.scale')
             n <- r$entryOrder
             palettes[[max(n)]][n]
           },
           'seq' = {
             ContinuousPtCol(seq_along(r$allTrees))
           },
           'dist' = {
             dat <- as.matrix(distances())[input$which.tree, ]
             message(paste0(round(dat, 2)[1:10], collapse = ', '))
             message(cut(dat, 256)[1:10])
             ContinuousPtCol(as.matrix(distances())[input$which.tree, ],
                             bigDark = TRUE)
           },
           'discrete' = {
             dat <- pointData()
             categories <- unique(dat)
             nCol <- length(categories)
             if (nCol <= length(palettes)) {
               PointDataStatus("")
               show('pt.col.scale')
               output$pt.col.scale <- renderPlot({
                 par(mar = c(1, 0, 0, 0))
                 plot(seq_len(nCol), rep(0, nCol), pch = 15, ylim = c(-1, 0),
                      col = palettes[[nCol]], ann = FALSE, axes = FALSE)
                 axis(1, at = seq_len(nCol), labels = categories,
                      tick = FALSE, line = -1)
               })
               palettes[[nCol]][match(dat, categories)]
             } else {
               hide('pt.col.scale')
               PointDataStatus(nCol, " categories is too many to display;",
                               "did you mean continuous?")
               c(palettes[[length(palettes)]], '#000000')[
                 pmin(match(dat, categories), length(palettes[[length(palettes)]]) + 1L)]
             }
           },
           'continuous' = {
             dat <- pointData()
             if (is.numeric(dat)) {
               ContinuousPtCol(dat)
             } else {
               hide('pt.col.scale')
               PointDataStatus("Point data are not numeric.",
                               'Try selecting "Custom categorical".')
               '#000000'
             }
           }
    )
  })
  
  pchs <- reactive({
    switch(input$pt.pch,
           '1' = 1,
           '16' = 16,
           'clust' = {
             cl <- clusterings()
             if (cl$sil > 0.25) {
               cl$cluster - 1
             } else 16
           },
           'entry' = {
             r$entryOrder - 1
           }
    )
  })
  
  ##############################################################################
  # Render plot
  ##############################################################################
  treespacePlot <- function() {
    cl <- clusterings()
    proj <- mapping()
    
    nDim <- min(dims(), nProjDim())
    plotSeq <- matrix(0, nDim, nDim)
    nPanels <- nDim * (nDim - 1L) / 2L
    plotSeq[upper.tri(plotSeq)] <- seq_len(nPanels)
    layout(t(plotSeq[-nDim, -1]))
    par(mar = rep(0.2, 4))
    withProgress(message = 'Drawing plot', {
      for (i in 2:nDim) for (j in seq_len(i - 1)) {
        incProgress(1 / nPanels)
        # Set up blank plot
        plot(proj[, j], proj[, i], ann = FALSE, axes = FALSE, frame.plot = TRUE,
             type = 'n', asp = 1, xlim = range(proj), ylim = range(proj))
        
        # Plot MST
        if (mstSize() > 0) {
          apply(mstEnds(), 1, function (segment)
            lines(proj[segment, j], proj[segment, i], col = "#bbbbbb", lty = 1))
        }
        
        # Add points
        points(proj[, j], proj[, i], pch = pchs(),
               col = paste0(pointCols(), as.hexmode(input$pt.opacity)),
               cex = input$pt.cex)
        
        if ("hulls" %in% input$display && cl$sil > 0.25) {
          # Mark clusters
          for (clI in seq_len(cl$n)) {
            inCluster <- cl$cluster == clI
            clusterX <- proj[inCluster, j]
            clusterY <- proj[inCluster, i]
            hull <- chull(clusterX, clusterY)
            polygon(clusterX[hull], clusterY[hull], lty = 1, lwd = 2,
                    border = palettes[[min(length(palettes), cl$n)]][clI])
            #border = '#54de25bb')
          }
        }
        if ('labelTrees' %in% input$display) {
          text(proj[, j], proj[, i], thinnedTrees())
        }
      }
    })
  }
  
  mode3D <- reactive("show3d" %in% input$display)
  PlotSize <- function () debounce(reactive(input$plotSize), 100)
  output$distPlot <- renderPlot({
    if (!mode3D()) {
      if (inherits(distances(), 'dist')) {
        treespacePlot()
        output$mappingStatus <- renderText('')
      } else {
        output$mappingStatus <- renderText("No distances available.")
      }
    }
  }, width = PlotSize(), height = PlotSize())
  
  output$threeDPlot <- rgl::renderRglwidget({
    if (mode3D() && inherits(distances(), 'dist')) {
      cl <- clusterings()
      proj <- mapping()
      withProgress(message = 'Drawing 3D plot', {
        rgl::rgl.open(useNULL = TRUE)
        incProgress(0.1)
        rgl::rgl.bg(color = 'white')
        rgl::plot3d(proj[, 1], proj[, 2], proj[, 3],
             aspect = 1, # Preserve aspect ratio - do not distort distances
             axes = FALSE, # Dimensions are meaningless
             col = pointCols(),
             alpha = input$pt.opacity / 255,
             cex = input$pt.cex,
             xlab = '', ylab = '', zlab = ''
        )
        incProgress(0.6)
        if ('labelTrees' %in% input$display) {
          rgl::text3d(proj[, 1], proj[, 2], proj[, 3], thinnedTrees())
        }
        if (mstSize() > 0) {
          apply(mstEnds(), 1, function (segment)
            rgl::lines3d(proj[segment, 1], proj[segment, 2], proj[segment, 3],
                    col = "#bbbbbb", lty = 1))
        }
      })
      rgl::rglwidget()
    }
  })
  
  output$savePng <- downloadHandler(
    filename = 'TreeSpace.png',
    content = function (file) {
      png(file, width = input$plotSize, height = input$plotSize)
      treespacePlot()
      dev.off()
  })
  
  output$savePdf <- downloadHandler(
    filename = 'TreeSpace.pdf',
    content = function (file) {
      pdf(file, title = paste0('Tree space mapping'))
      treespacePlot()
      dev.off()
  })
  
  output$saveCons <- downloadHandler(
    filename = 'ClusterConsensusTrees.pdf',
    content = function (file) {
      cl <- clusterings()
      pdf(file, title = if (cl$sil > 0.25) {
        paste0("Consensus trees for ", cl$n, " clusters found with ", cl$method,
               " using ", chosenDistance())
      } else {
        paste0("Consensus of all trees (no meaningful clusters found using ",
               chosenDistance(), ")")
      },
      width = 8, height = consSize() / 100, pointsize = 10)
      PlotClusterCons()
      dev.off()
  })
  
  
  
  ##############################################################################
  # References
  ##############################################################################
  output$references <- renderUI({
    tags$div(style = paste0('position: relative; top: ', 
                            (input$plotSize - 600)
                            + if ('cons' %in% input$display) consSize() - 200 else 0
                            , 'px'),
    tagList(
      tags$h2('References for methods used'),
      #HTML(paste0("<h2>References", input$plotSize, "</h2>")),
      tags$h3('Tree space construction'),
      HTML(paste0(Smith2021,
                  Kaski2003, Venna2001, RCoreTeam)),
      HTML(if (mstSize() > 0) Gower1969),
      HTML(if(input$distance == 'qd') SmithQuartet),
      tags$h3('Tree distance'),
      HTML(switch(input$distance,
             'cid' = paste0(Smith2020, SmithDist),
             'pid' = paste0(Smith2020, SmithDist),
             'qd' = paste0(Estabrook1985, Sand2014, SmithQuartet),
             'kc' = paste0(Kendall2016),
             'path' = paste0(Farris1973, Steel1993, SmithDist),
             'rf' = paste0(Robinson1981, Day1985, SmithDist))
      ),
      tags$h3('Mapping'),
      HTML(switch(input$mapping,
                  'pca' = Gower1966,
                  'k' = paste0(Kruskal1964, Venables2002),
                  'nls' = paste0(Sammon1969, Venables2002)
                  )
           ),
      tags$h3('Clustering'),
      HTML(paste(c(pam = paste0('Partitioning around medoids: ', Maechler2019),
        hmm = paste0("Hierarchical, minimax linkage: ", Bien2011, Murtagh1983),
        hsi = '',#paste0("Hierarchical, single linkage:"),
        hco = '',#paste0("Hierarchical, complete linkage:"),
        hav = '',#paste0("Hierarchical, average linkage:"),
        hmd = '',#paste0("Hierarchical, median linkage:"),
        hct = '',#paste0("Hierarchical, centroid linkage:"),
        hwd = paste0("Hierarchical, Ward d\ub2 linkage: ", Ward1963),
        kmn = '',#paste0("K-means:"),
        spec = ''#paste0("Spectral:")
        )[input$clustering])),
      HTML(paste("Cluster consensus trees:", Stockham2002))
    ))
  })
}

shinyApp(ui = ui, server = server)
