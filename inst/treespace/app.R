library("shiny")
library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)

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
                 
                 

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Tree space analysis"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      fileInput("treefile", "Trees", placeholder = "No tree file selected"),
      textOutput(outputId = "treesStatus"),
      textOutput(outputId = "keptTrees"),
      sliderInput(inputId = "keepTrees",
                  label = "Crop ends (%):",
                  min = 0,
                  max = 100,
                  value = c(0, 100)),
      textOutput(outputId = "thinnedTrees"),
      sliderInput(inputId = "thinTrees",
                  label = "Sample every 2^:",
                  min = 0,
                  max = 12,
                  round = FALSE,
                  step = 0.01,
                  value = 1),
      actionButton("addTrees", "Add to existing"),
      actionButton("replaceTrees", "Replace existing"),
      textOutput(outputId = "treesInMemory"),
                                        
      radioButtons("distance", "Distance method",
                   choices = list("Clustering information" = 'cid',
                                  "Phylogenetic information" = 'pid',
                                  "Quartet" = 'qd',
                                  "Path" = 'path',
                                  "Robinson\u2212Foulds" = 'rf'),
                   selected = 'cid'),
      textOutput(outputId = "distanceStatus"),
      
      
      radioButtons("projection", "Projection method",
                   choices = list("Principal components (classical MDS)" = 'pca',
                                  "Kruskal-1 nmMDS" = 'k',
                                  "Sammon mapping mMDS" = 'nls'),
                   selected = 'pca'),
      textOutput(outputId = "projectionStatus"),
      checkboxInput('labelTrees', 'Plot tree numbers', FALSE),
      
      sliderInput(inputId = "dims",
                  label = "Dimensions to plot:",
                  min = 1,
                  max = 15,
                  value = 6),
      
      sliderInput(inputId = "mst",
                  label = "Minimum spanning tree extent:",
                  min = 0,
                  max = 100,
                  value = 20),
      
      checkboxGroupInput("clustering", "Clustering methods",
                         choices = list("Partitioning around medoids" = 'pam',
                                        "Heirarchical, minimax linkage" = 'hmm',
                                        "Heirarchical, TBC linkage" = 'har',
                                        "K-means" = 'k',
                                        "Spectral" = "spec"),
                         selected = c('pam', 'hmm')),
      sliderInput(inputId = "maxClusters",
                  label = "Maximum cluster number:",
                  min = 2,
                  max = 20,
                  value = 8),
      textOutput(outputId = "clusteringStatus"),
      
    ),
    
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: plot
      plotOutput(outputId = "distPlot", height = "800px"),
      plotOutput(outputId = "projPlot",  width = "45%", height = "300px"),
      plotOutput(outputId = "clustPlot", width = "45%", height = "300px"),
      
    )
  ),
  
  # fluidRow(
  # )
  
)

server <- function(input, output) {
  
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
    if (is.null(input)) return("No tree file selected")
    tmpFile <- fileInput$datapath
    if (is.null(tmpFile)) return ("No trees found")
    if (length(grep('#NEXUS', toupper(readLines(tmpFile)[1]),
                    fixed = TRUE)) > 0) {
      ret <- read.nexus(tmpFile)
    } else {
      ret <- read.tree(tmpFile)
    }
    
    if (!inherits(ret, c('phylo', 'multiPhylo'))) {
      return("Could not read trees from file")
    }
    ret
  })

  output$treesStatus <- renderText({
    if (is.character(treeFile())) return(treeFile())
    nTrees <- length(treeFile())
    paste(nTrees, "trees in file.")
  })
  
  keptRange <- reactive({
    nTrees <- length(treeFile())
    keeps <- input$keepTrees
    c(max(1L, ceiling(keeps[1] * nTrees / 100)),
      floor(keeps[2] * nTrees / 100))
  })
  
  output$keptTrees <- renderText({
    nTrees <- length(treeFile())
    paste0("Keeping trees ", keptRange()[1]," to ", keptRange()[2], '.')
  })
  
  output$thinnedTrees <- renderText({
    if (input$thinTrees == 1) "" else {
      paste("Sampled", length(thinnedTrees()), "trees.")
    }
  })

  thinnedTrees <- reactive({
    seq(keptRange()[1], keptRange()[2], by = 2^input$thinTrees)
  })
  
  observeEvent(input$addTrees, {
    newTrees <- treeFile()[thinnedTrees()]
    r$allTrees <- if (length(r$allTrees) == 0) newTrees else {
      c(r$allTrees, newTrees)
    }
    ClearDistances()
  })
  
  observeEvent(input$replaceTrees, {
    r$allTrees <- c(treeFile()[thinnedTrees()])
    ClearDistances()
  })
  
  output$treesInMemory <- renderText({
    paste(length(r$allTrees), "trees in memory.")
  })
  
  distances <- reactive({
    switch(input$distance,
           'cid' = {
             if (is.null(r$dist_cid)) r$dist_cid = ClusteringInfoDistance(r$allTrees)
             output$distanceStatus <- renderText('CI distances calculated.')
             r$dist_cid
           },
           'pid' = {
             if (is.null(r$dist_pid)) r$dist_pid = PhylogeneticInfoDistance(r$allTrees)
             output$distanceStatus <- renderText('Phylo. Info. distances calculated.')
             r$dist_pid
           },
           'qd' = {
             if (is.null(r$dist_qd)) r$dist_qd = as.dist(Quartet::QuartetDivergence(
               Quartet::ManyToManyQuartetAgreement(r$allTrees), similarity = FALSE))
             output$distanceStatus <- renderText('Quartet distances calculated.')
             r$dist_qd
           }, 'path' = {
             if (is.null(r$dist_path)) r$dist_path = PathDist(r$allTrees)
             output$distanceStatus <- renderText('Path distances calculated.')
             r$dist_path
           }, 'rf' = {
             if (is.null(r$dist_rf)) r$dist_rf = RobinsonFoulds(r$allTrees)
             output$distanceStatus <- renderText('RF distances calculated.')
             r$dist_rf
           }
    )
    
  })
  
  projection <- reactive({
    proj_id <- paste0('proj_', input$projection, '_', input$distance)
    if (is.null(r[[proj_id]])) {
      proj <- switch(input$projection,
                     'pca' = cmdscale(distances(), k = 15),
                     'k' = MASS::isoMDS(distances(), k = 15)$points,
                     'nls' = MASS::sammon(distances(), k = 15)$points
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
      
      hSil <- pamSil <- harSil <- kSil <- specSil <- -99
      dists <- distances()
      
      if ('pam' %in% input$clustering) {
        pamClusters <- lapply(possibleClusters, function (k) {
          cluster::pam(dists, k = k)
        })
        pamSils <- vapply(pamClusters, function (pamCluster) {
          mean(cluster::silhouette(pamCluster)[, 3])
        }, double(1))
      
        bestPam <- which.max(pamSils)
        pamSil <- pamSils[bestPam]
        pamCluster <- pamClusters[[bestPam]]$cluster
      }
      
      if ('hmm' %in% input$clustering) {
        hTree <- protoclust::protoclust(dists)
        hClusters <- lapply(possibleClusters, function (k) cutree(hTree, k = k))
        hSils <- vapply(hClusters, function (hCluster) {
          mean(cluster::silhouette(hCluster, dists)[, 3])
        }, double(1))
        bestH <- which.max(hSils)
        hSil <- hSils[bestH]
        hCluster <- hClusters[[bestH]]
      }
      bestCluster <- c('pam', 'hmm', 'har', 'k', 'spec')[
        which.max(c(pamSil, hSil, harSil, kSil, specSil))]
      r[[clust_id]] <- list(method = switch(bestCluster,
                                            pam = 'partitioning around medoids',
                                            hmm = 'minimax linkage',
                                            har = 'average linkage',
                                            k = 'K-means',
                                            spec = 'spectral'),
                            n = 1 + switch(bestCluster,
                                           pam = bestPam, hmm = bestH, 
                                           har = bestHav, k = bestK, spec = bestSpec),
                            sil = switch(bestCluster,
                                         pam = pamSil, hmm = hSil,
                                         har = havSil, k = kSil, spec = specSil), 
                            cluster = switch(bestCluster,
                                             pam = pamCluster,
                                             hmm = hCluster,
                                             har = havCluster,
                                             k = kCluster,
                                             spec = specCluster)
      )
    }
    r$clust_max <- maxCluster
    r[[clust_id]]
  })
  
  output$clusteringStatus <- renderText({
    cl <- clusterings()
    paste0(cl$n, " clusters found with ", cl$method, '; silhouette coeff. = ', 
           round(cl$sil, 3))
  })
  
  mstEnds <- reactive({
    
    mstSize <- input$mst
    if (!is.null(r$mst_size) && mstSize != r$mst_size) ClearMST()
    mst_id <- paste0('mst_', input$distance)
    
    if (is.null(r[[mst_id]])) {
      dist <- as.matrix(distances())
      nRows <- dim(dist)[1]
      selection <- unique(round(seq.int(1, nRows, length.out = max(2L, mstSize / 100 * nRows))))
      edges <- MSTEdges(dist[selection, selection])
      edges[, 1] <- selection[edges[, 1]]
      edges[, 2] <- selection[edges[, 2]]
      r[[mst_id]] <- edges
    }
    
    r$mst_size <- mstSize
    r[[mst_id]]
  })
  
  output$distPlot <- renderPlot({
    output$projectionStatus <- renderText('Projecting...')
    
    if (inherits(distances(), 'dist')) {
      cl <- clusterings()
      col <- if (cl$sil > 0.25) {
        palettes[[min(length(palettes), cl$n)]][cl$cluster]
      } else palettes[[1]]
      proj <- projection()
      plot(proj,
           asp = 1, # Preserve aspect ratio - do not distort distances
           ann = FALSE, axes = FALSE, # Dimensions are meaningless
           frame.plot = TRUE,
           pch = 16,
           col = col
      )
      if (input$labelTrees) text(proj[, 1], proj[, 2], seq_len(nrow(proj)))
      
      if (input$mst > 0) {
        apply(mstEnds(), 1, function (segment)
          lines(proj[segment, 1], proj[segment, 2], col = "#bbbbbb", lty = 1))
      }
      
      output$projectionStatus <- renderText(paste(switch(input$projection, 
                                                   'pca' = 'Principal components',
                                                   'k' = 'Kruskal-1',
                                                   'nls' = 'Sammon'),
                                            'projection plotted.'))
    } else {
      output$projectionStatus <- renderText("No distances available.")
    }
  })
  
  output$clustPlot <- renderPlot({
    
    #plot(pamSils ~ possibleClusters,
    #     xlab = 'Number of clusters', ylab = 'Silhouette coefficient',
    #     ylim = range(c(pamSils, hSils)))
    #points(hSils ~ possibleClusters, pch = 2)
    #legend('topright', c('PAM', 'Hierarchical'), pch = 1:2)
  })
  
}

shinyApp(ui = ui, server = server)