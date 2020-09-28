library("shiny")
library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
library("TreeDist")

options(shiny.maxRequestSize = 100 * 1024^2)

treeNumbers <- c(1:220)
trees <- as.phylo(treeNumbers, 8)
spectrum <- viridisLite::plasma(220)
treeCols <- spectrum[treeNumbers]

palettes <- list("#91aaa7",
                 c("#969660", "#c3dfca"),
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

# Define UI for app that draws a histogram ----
ui <- fluidPage(theme = 'treespace.css',
  
  # Sidebar layout with input and output definitions ----
  
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
                      value = c(0, 90)),
          sliderInput(inputId = "thinTrees",
                      label = "Sample every:",
                      min = 0,
                      max = 12,
                      round = FALSE,
                      step = 0.01,
                      pre = '2^',
                      value = 1),
          textOutput(outputId = "thinnedTrees"),
          actionButton("replaceTrees", "Replace existing"),
          actionButton("addTrees", "Add to existing"),
        ),
        tabPanel('Analysis',
           selectInput("distance", "Distance method",
                       choices = list("Clustering information" = 'cid',
                                      "Phylogenetic information" = 'pid',
                                      "Quartet (slow)" = 'qd',
                                      "Path" = 'path',
                                      "Robinson\u2212Foulds" = 'rf'),
                       selected = 'cid'),
          
          textOutput(outputId = "projectionStatus"),
          fluidRow(plotOutput(outputId = "pcQuality", height = "72px")),
          selectInput("projection", "Projection method",
                       choices = list("Princ. comps (classical MDS)" = 'pca',
                                      "Sammon mapping mMDS" = 'nls',
                                      "Kruskal-1 nmMDS (slow)" = 'k'
                                      ),
                       selected = 'pca'),
          
          checkboxGroupInput("clustering", "Clustering methods",
                             choices = list("Partitioning around medoids" = 'pam',
                                            "Heirarchical, minimax linkage" = 'hmm',
                                            "Heirarchical, single linkage" = 'hsi',
                                            "Heirarchical, complete linkage" = 'hco',
                                            "Heirarchical, average linkage" = 'hav',
                                            "Heirarchical, median linkage" = 'hmd',
                                            "Heirarchical, centroid linkage" = 'hct',
                                            "Heirarchical, Ward d\ub2 linkage" = 'hwd',
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
        )
      ),
    ),
    
    column(9,
      mainPanel(
        fluidRow(id = 'plotConfig',
          tags$span("Plot size:", id = 'plotSizeSpan'),
          sliderInput(inputId = "plotSize", label = NULL,
                      width = '200px',
                      min = 100, max = 2000,
                      post = 'px', value = 600),
          tags$span("Save as: "),
          downloadButton('savePdf', 'PDF'),
          downloadButton('savePng', 'PNG'),
        ),
        fluidRow(
          plotOutput(outputId = "distPlot"),
          rgl::rglwidgetOutput(outputId = "threeDPlot", width = "300px", height = "300px"),
          ),
      )
  )
  
)

server <- function(input, output, session) {
  
  #r <- reactiveValues(allTrees = structure(list(), class = 'multiPhylo'))
  #treeNumbers <- c(1:220)
  treeNumbers <- c(1:50, 401:440)
  
  r <- reactiveValues(allTrees = as.phylo(treeNumbers, 8))
  
  ClearClusters <- function () {
    r$clust_cid = NULL
    r$clust_pid = NULL
    r$clust_qd = NULL
    r$clust_path = NULL
    r$clust_rf = NULL
  }
  
  ClearMST <- function () {
    r$mst_cid = NULL
    r$mst_pid = NULL
    r$mst_qd = NULL
    r$mst_path = NULL
    r$mst_rf = NULL
  }
  
  ClearDistances <- function () {
    r$dist_cid = NULL
    r$dist_pid = NULL
    r$dist_qd = NULL
    r$dist_path = NULL
    r$dist_rf = NULL
    
    r$proj_pca_cid = NULL
    r$proj_pca_pid = NULL
    r$proj_pca_qd = NULL
    r$proj_pca_path = NULL
    r$proj_pca_rf = NULL
    
    r$proj_k_cid = NULL
    r$proj_k_pid = NULL
    r$proj_k_qd = NULL
    r$proj_k_path = NULL
    r$proj_k_rf = NULL
    
    r$proj_nls_cid = NULL
    r$proj_nls_pid = NULL
    r$proj_nls_qd = NULL
    r$proj_nls_path = NULL
    r$proj_nls_rf = NULL
    
    ClearClusters()
    
    ClearMST()
  }
  
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
  
  observeEvent(input$addTrees, {
    newTrees <- treeFile()[thinnedTrees()]
    r$allTrees <- if (length(r$allTrees) == 0) newTrees else {
      c(r$allTrees, newTrees)
    }
    nTrees <- length(r$allTrees)
    updateSliderInput(session, 'mst', value = min(100, nTrees), max = nTrees)
    ClearDistances()
  })
  
  observeEvent(input$replaceTrees, {
    r$allTrees <- c(treeFile()[thinnedTrees()])
    nTrees <- length(r$allTrees)
    updateSliderInput(session, 'mst', value = min(100, nTrees), max = nTrees)
    ClearDistances()
  })
  
  output$treesInMemory <- renderText({
    paste(length(r$allTrees), "trees in memory.")
  })
  
  distances <- reactive({
    withProgress(message = 'Calculating distances', value = 0.99,
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
           }, 'path' = {
             if (is.null(r$dist_path)) r$dist_path = PathDist(r$allTrees)
             r$dist_path
           }, 'rf' = {
             if (is.null(r$dist_rf)) r$dist_rf = RobinsonFoulds(r$allTrees)
             r$dist_rf
           }
    ))
    
  })
  
  projection <- reactive({
    proj_id <- paste0('proj_', input$projection, '_', input$distance)
    if (is.null(r[[proj_id]])) {
      withProgress(
        message = 'Projecting distances',
        value = 0.99,
        proj <- switch(input$projection,
                       'pca' = cmdscale(distances(), k = 15),
                       'k' = MASS::isoMDS(distances(), k = 15)$points,
                       'nls' = MASS::sammon(distances(), k = 15)$points
                       )
      )
      r[[proj_id]] <- proj
    }
    r[[proj_id]]
  })
  
  
  clusterings <- reactive({
    maxCluster <- input$maxClusters
    if (!is.null(r$clust_max) && maxCluster != r$clust_max) ClearClusters()
    clust_id <- paste0('clust_', input$distance)
    
    if (is.null(r[[clust_id]])) {
      possibleClusters <- 2:input$maxClusters
      
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
          hTree <- stats::hclust(dists, method = 'centroid')
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
          spectralEigens <- SpectralClustering(dists,
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
      })
    }
    r[[clust_id]]
  })
  
  observeEvent(input$calcClusters, {ClearClusters()})
  
  mstEnds <- reactive({
    
    mstSize <- input$mst
    if (!is.null(r$mst_size) && mstSize != r$mst_size) ClearMST()
    mst_id <- paste0('mst_', input$distance)
    if (is.null(r[[mst_id]])) {
      dist <- as.matrix(distances())
      nRows <- dim(dist)[1]
      withProgress(message = 'Calculating MST', {
        selection <- unique(round(seq.int(1, nRows, length.out = max(2L, mstSize))))
        edges <- MSTEdges(dist[selection, selection])
        edges[, 1] <- selection[edges[, 1]]
        edges[, 2] <- selection[edges[, 2]]
        r[[mst_id]] <- edges
      })
    }
    
    r$mst_size <- mstSize
    r[[mst_id]]
  })
  
  treespacePlot <- function () {
    cl <- clusterings()
    treeCols <- if (cl$sil > 0.25) {
      palettes[[min(length(palettes), cl$n)]][cl$cluster]
    } else palettes[[1]]
    proj <- projection()
    
    nDim <- input$dims
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
        if (input$mst > 0) {
          apply(mstEnds(), 1, function (segment)
            lines(proj[segment, j], proj[segment, i], col = "#bbbbbb", lty = 1))
        }
        
        # Add points
        points(proj[, j], proj[, i], pch = 16,
               col = paste0(treeCols, as.hexmode(input$pt.opacity)),
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

  PlotSize <- function () debounce(reactive(input$plotSize), 100)
  output$distPlot <- renderPlot({
    if (inherits(distances(), 'dist')) {
      treespacePlot()
      output$projectionStatus <- renderText('')
    } else {
      output$projectionStatus <- renderText("No distances available.")
    }
  }, width = PlotSize(), height = PlotSize())
  
  output$threeDPlot <- rgl::renderRglwidget({
    if ("show3d" %in% input$display) {
      
      if (inherits(distances(), 'dist')) {
        cl <- clusterings()
        col <- if (cl$sil > 0.25) {
          palettes[[min(length(palettes), cl$n)]][cl$cluster]
        } else palettes[[1]]
        proj <- projection()
        withProgress(message = 'Drawing 3D plot', {
          rgl::rgl.open(useNULL = TRUE)
          incProgress(0.1)
          rgl::rgl.bg(color = 'white')
          rgl::plot3d(proj[, 1], proj[, 2], proj[, 3],
               aspect = 1, # Preserve aspect ratio - do not distort distances
               axes = FALSE, # Dimensions are meaningless
               pch = 16,
               col = treeCols,
               alpha = input$pt.opacity / 255,
               cex = input$pt.cex,
               xlab = '', ylab = '', zlab = ''
          )
          incProgress(0.6)
          if ('labelTrees' %in% input$display) {
            rgl::text3d(proj[, 1], proj[, 2], proj[, 3], thinnedTrees())
          }
          if (input$mst > 0) {
            apply(mstEnds(), 1, function (segment)
              rgl::lines3d(proj[segment, 1], proj[segment, 2], proj[segment, 3],
                      col = "#bbbbbb", lty = 1))
          }
        })
        rgl::rglwidget()
      }
    }
  })
  
  dims <- debounce(reactive(input$dims), 400)
  projQual <- reactive({
    withProgress(message = "Estimating projection quality", {
      vapply(1:15, function (k) {
        incProgress(1/15)
        ProjectionQuality(distances(), dist(projection()[, seq_len(k)]), 10)
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
    
    logScore <- LogScore(projQual()['TxC', input$dims])
    lines(rep(logScore, 2), c(-1, 1), lty = 3)
    
    tickPos <- c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0)
    ticks <- LogScore(tickPos)
    
    axis(1, at = ticks, labels = NA, line = 0)
    axis(1, tick = FALSE, at = ticks, labels = tickPos, line = 0)
    axis(1, line = -1, tick = FALSE,
         at = ticks[-1] - ((ticks[-1] - ticks[-length(ticks)]) / 2),
         labels = c('', 'dire', '', "ok", "gd", "excellent"))
    axis(3, at = 0.5, tick = FALSE, line = -2, 
         paste0(dims(), 'D projection quality (trustw. \ud7 contin.):'))
  })
  
  output$clustQuality <- renderPlot({
    par(mar = c(2, 0.5, 0, 0.5), xpd = NA, mgp = c(2, 1, 0))
    cl <- clusterings()
    clust_id <- paste0('clust_', input$distance)
    sil <- r[[clust_id]]$sil
    if (length(sil) == 0) sil <- -0.5
    nStop <- 400
    range <- c(0.5, 1)
    negScale <- viridisLite::plasma(nStop)[seq(range[1] * nStop, 1,
                                               length.out = nStop * range[1])]
    posScale <- viridisLite::viridis(nStop)
    
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
         labels = paste0("Silhouette coefficient (", round(cl$sil, 3), ")"))
         
    
    axis(3, tick = FALSE, line = -2, at = 0.25,
         labels = paste0(cl$n, " clusters found with ", cl$method))
         
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
      pdf(file, title = paste0('Tree space projection'))
      treespacePlot()
      dev.off()
  })
}

shinyApp(ui = ui, server = server)