library("shiny")
library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)

treeNumbers <- c(1:220)
trees <- as.phylo(treeNumbers, 8)
spectrum <- viridisLite::plasma(220)
treeCols <- spectrum[treeNumbers]

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
      
      checkboxGroupInput("clustering", "Clustering methods",
                   choices = list("Partitioning around medoids" = 'pam',
                                  "Heirarchical, minimax linkage" = 'hmm',
                                  "Heirarchical, TBC linkage" = 'har',
                                  "K-means" = 'k',
                                  "Spectral" = "spec"),
                   selected = c('pam', 'hmm')),
      actionButton("calcClusters", "Calculate clusterings"),
      textOutput(outputId = "clusteringStatus"),
      
      radioButtons("projection", "Projection method",
                   choices = list("Principal components (classical MDS)" = 'pca',
                                  "Kruskal-1 nmMDS" = 'k',
                                  "Sammon mapping mMDS" = 'nls'),
                   selected = 'pca'),
      textOutput(outputId = "projectionStatus"),
      
      plotOutput(outputId = "projectionStress"),
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "dims",
                  label = "Dimensions to plot:",
                  min = 1,
                  max = 15,
                  value = 6)
      
    ),
    
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: plot
      plotOutput(outputId = "distPlot", height = "800px")
      
    )
  ),
  
  # fluidRow(
  # )
  
)

server <- function(input, output) {
  
  #r <- reactiveValues(allTrees = structure(list(), class = 'multiPhylo'))
  treeNumbers <- c(1:220)
  
  r <- reactiveValues(allTrees = as.phylo(treeNumbers, 8))
  
  ClearDistances <- function () {
    r$dist_cid = NULL
    r$dist_pid = NULL
    r$dist_qd = NULL
    r$dist_path = NULL
    r$dist_rf = NULL
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
  
  
  output$distPlot <- renderPlot({
    output$projectionStatus <- renderText('Projecting...')
    
    if (inherits(distances(), 'dist')) {
      projection <- switch(input$projection, 
             'pca' = cmdscale(distances(), k = 15),
             'k' = MASS::isoMDS(distances(), k = 15)$points,
             'nls' = MASS::sammon(distances(), k = 15)$points)
      plot(projection,
           asp = 1, # Preserve aspect ratio - do not distort distances
           ann = FALSE, axes = FALSE, # Dimensions are meaningless
           frame.plot = TRUE,
           pch = 16
      )
      output$projectionStatus <- renderText(paste(switch(input$projection, 
                                                   'pca' = 'Principal components',
                                                   'k' = 'Kruskal-1',
                                                   'nls' = 'Sammon'),
                                            'projection'))
    } else {
      output$projectionStatus <- renderText("No distances available.")
    }
  })
  
}

shinyApp(ui = ui, server = server)