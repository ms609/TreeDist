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
      sliderInput(inputId = "keepTrees",
                  label = "Crop ends (%):",
                  min = 0,
                  max = 100,
                  value = c(0, 100)),
      sliderInput(inputId = "thinTrees",
                  label = "Sample every 2^:",
                  min = 0,
                  max = 12,
                  round = FALSE,
                  step = 0.01,
                  value = 1),
      actionButton("addTrees", "Add to existing"),
      actionButton("replaceTrees", "Replace existing"),
                                        
      checkboxGroupInput("distances", "Distance methods",
                         choices = list("clustering information" = 'cid',
                                        "phylogenetic information" = 'pid',
                                        "quartet" = 'qd',
                                        "path" = 'path',
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
      plotOutput(outputId = "distPlot")
      
    )
  ),
  
  # fluidRow(
  # )
  
)

server <- function(input, output) {
  
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
  
  allTrees <- reactive({
    
    #if (inherits(tmp, 'multiPhylo'))
    if (input$addTrees) {
      c(trees, tmp)
    } else {
      tmp
    }
    class(tmp)
  })
  
  keptRange <- reactive({
    nTrees <- length(treeFile())
    keeps <- input$keepTrees
    c(max(1L, ceiling(keeps[1] * nTrees / 100)),
      floor(keeps[2] * nTrees / 100))
  })
  
  thinnedTrees <- reactive({
    seq(keptRange()[1], keptRange()[2], by = 2^input$thinTrees)
  })
  
  output$treesStatus <- renderText({
    if (is.character(treeFile())) return(treeFile())
    nTrees <- length(treeFile())
    found <- paste("Found", nTrees, "trees in file.")
    keeps <- if (keptRange()[1] == 1 && keptRange()[2] == nTrees) "" else {
      paste0("Keeping trees ", keptRange()[1]," to ", keptRange()[2], '.')
    }
    sampling <- if (input$thinTrees == 1) "" else {
      paste("Uniformly sampling", length(thinnedTrees()), "trees.")
    }
    
    paste(found, keeps, sampling)
  })
  
  distances_rf <- reactive({RobinsonFoulds(trees)})
  
  
  output$distPlot <- renderPlot({
    
    
  })
  
}

shinyApp(ui = ui, server = server)