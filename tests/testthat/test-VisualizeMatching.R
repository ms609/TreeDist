library("TreeTools")

test_that("#86 is fixed", {
  tree1 <- read.tree(text = "((((GCF_002862005_1_ASM286200v1_genomic:0.0,GCF_002861945_1_ASM286194v1_genomic:0.0,GCF_000213955_1_ASM21395v1_genomic:0.0,GCF_013315085_1_ASM1331508v1_genomic:0.0,GCF_002861975_1_ASM286197v1_genomic:0.0,GCF_013315025_1_ASM1331502v1_genomic:0.0,GCF_002861965_1_ASM286196v1_genomic:0.0,GCF_002862015_1_ASM286201v1_genomic:0.0,GCF_013315045_1_ASM1331504v1_genomic:0.0):0.000000006,((GCF_000414525_1_ASM41452v1_genomic:0.002230831,(((GCF_001546445_1_ASM154644v1_genomic:0.0,GCF_013315115_1_ASM1331511v1_genomic:0.0):0.000000005,((GCF_002861905_1_ASM286190v1_genomic:0.0,GCF_000414605_1_ASM41460v1_genomic:0.0,GCF_000414665_1_ASM41466v1_genomic:0.0,GCF_000414585_1_ASM41458v1_genomic:0.0):0.000000005,GCF_002861885_1_ASM286188v1_genomic:0.002230870)0.966:0.006820446)0.969:0.009438344,((GCF_003426565_1_ASM342656v1_genomic:0.0,piotii_GCF_003397585_1_ASM339758v1_genomic:0.0):0.000000005,(((GCF_000414545_1_ASM41454v1_genomic:0.003992352,(GCF_003408835_1_ASM340883v1_genomic:0.0,GCF_000414505_1_ASM41450v1_genomic:0.0):0.000000005)0.909:0.008040404,(GCF_003397615_1_ASM339761v1_genomic:0.000000005,(GCF_000414565_1_ASM41456v1_genomic:0.003708023,(GCF_000414625_1_ASM41462v1_genomic:0.016167672,(GCF_000414485_1_ASM41448v1_genomic:0.0,GCF_000414425_1_ASM41442v1_genomic:0.0):0.000000005)0.000:0.000000005)0.000:0.000000006)0.928:0.000000005)0.948:0.018701881,(GCF_001546455_1_ASM154645v1_genomic:0.048525006,(GCF_001563665_1_ASM156366v1_genomic:0.116224039,((GCF_003408775_1_ASM340877v1_genomic:0.000000005,GCF_002884775_1_ASM288477v1_genomic:0.015930854)0.913:0.012725504,(((GCF_001953155_1_ASM195315v1_genomic:0.001988946,(GCF_013315145_1_ASM1331514v1_genomic:0.000000005,((GCF_003397745_1_ASM339774v1_genomic:0.0,GCF_003408815_1_ASM340881v1_genomic:0.0,swidsinskii_GCF_003397705_1_ASM339770v1_genomic:0.0):0.000000005,(GCF_000025205_1_ASM2520v1_genomic:0.000000005,(GCF_002884855_1_ASM288485v1_genomic:0.0,GCF_002884875_1_ASM288487v1_genomic:0.0):0.001978715)0.931:0.003992623)0.000:0.000000005)0.469:0.000000006)0.885:0.005473514,(GCF_013315195_1_ASM1331519v1_genomic:0.002012914,((GCF_002861125_1_ASM286112v1_genomic:0.0,GCF_013315125_1_ASM1331512v1_genomic:0.0,GCF_013315255_1_ASM1331525v1_genomic:0.0,GCF_003397635_1_ASM339763v1_genomic:0.0,GCF_002861145_1_ASM286114v1_genomic:0.0):0.006168310,leopoldii_GCF_003293675_1_ASM329367v1_genomic:0.001998855)0.781:0.002009527)0.871:0.004589922)0.934:0.014500618,(GCF_003408845_1_ASM340884v1_genomic:0.033674519,((GCF_000414465_1_ASM41446v1_genomic:0.0,GCF_000414445_1_ASM41444v1_genomic:0.0):0.001943013,GCF_001546485_1_ASM154648v1_genomic:0.000000005)1.000:0.047170589)0.714:0.015802120)0.924:0.021101316)0.654:0.021430364)0.995:0.072887040)0.435:0.003580847)0.278:0.005651482)0.892:0.005483145)0.884:0.005371297)0.793:0.000000005,GCF_003408785_1_ASM340878v1_genomic:0.004504993)0.849:0.002683606)0.000:0.000000005,(GCF_001660735_1_ASM166073v1_genomic:0.0,GCF_013315075_1_ASM1331507v1_genomic:0.0):0.000000005)0.932:0.005409354,(GCF_002861165_1_ASM286116v1_genomic:0.0,GCF_001660755_1_ASM166075v1_genomic:0.0):0.001866731,((GCF_000414645_1_ASM41464v1_genomic:0.000000005,(GCF_900637625_1_52295_C01_genomic:0.0,GCF_000414685_1_ASM41468v1_genomic:0.0,GCF_001042655_1_ASM104265v1_genomic:0.0,GCF_003397665_1_ASM339766v1_genomic:0.0,GCF_003408745_1_ASM340874v1_genomic:0.0,GCF_000159155_2_ASM15915v2_genomic:0.0,GCF_000178355_1_ASM17835v1_genomic:0.0,GCF_013315005_1_ASM1331500v1_genomic:0.0):0.000000005)0.489:0.000000005,((GCF_003585655_1_ASM358565v1_genomic:0.0,GCF_000414705_1_ASM41470v1_genomic:0.0,GCF_003812765_1_ASM381276v1_genomic:0.0,GCF_003585755_1_ASM358575v1_genomic:0.0):0.000000005,GCF_003397605_1_ASM339760v1_genomic:0.004518743)0.000:0.000000005)0.748:0.000000005);")
  text2 <-  paste0(
    "(GCF_013315005_1_ASM1331500v1_genomic:0.029801777,((GCF_002861945_1_ASM286194v1_genomic:0.000000005,GCF_013315085_1_ASM1331508v1_genomic:0.000031313)1.000:0.021319502,(GCF_000159155_2_ASM15915v2_genomic:0.000031236,(GCF_000178355_1_ASM17835v1_genomic:0.000438654,(GCF_001042655_1_ASM104265v1_genomic:0.000062623,GCF_900637625_1_52295_C01_genomic:0.000031311)0.387:0.000000005)0.928:0.000094032)1.000:0.021525169)1.000:0.008877453,((((GCF_003585655_1_ASM358565v1_genomic:0.037400199,(((GCF_001563665_1_ASM156366v1_genomic:0.652791246,(((GCF_002884775_1_ASM288477v1_genomic:0.079035837,GCF_003408775_1_ASM340877v1_genomic:0.090527862)1.000:0.099122682,((leopoldii_GCF_003293675_1_ASM329367v1_genomic:0.010984029,(GCF_003397635_1_ASM339763v1_genomic:0.010125736,((GCF_002861125_1_ASM286112v1_genomic:0.000000005,GCF_002861145_1_ASM286114v1_genomic:0.000031301)1.000:0.009693544,((GCF_013315125_1_ASM1331512v1_genomic:0.0,GCF_013315255_1_ASM1331525v1_genomic:0.0):0.012966743,GCF_013315195_1_ASM1331519v1_genomic:0.013648178)1.000:0.005988651)1.000:0.004698554)0.990:0.005500145)1.000:0.069052794,((GCF_002884855_1_ASM288485v1_genomic:0.000062668,GCF_002884875_1_ASM288487v1_genomic:0.000000005)1.000:0.025579950,(((GCF_013315145_1_ASM1331514v1_genomic:0.019646990",
    ",swidsinskii_GCF_003397705_1_ASM339770v1_genomic:0.029512978)1.000:0.013756732,GCF_003408815_1_ASM340881v1_genomic:0.023520587)0.995:0.005875475,(GCF_003397745_1_ASM339774v1_genomic:0.030659080,(GCF_000025205_1_ASM2520v1_genomic:0.025213714,GCF_001953155_1_ASM195315v1_genomic:0.028312469)1.000:0.011831021)0.833:0.005198931)1.000:0.031300458)1.000:0.042156286)1.000:0.115485782)1.000:0.022843945,((GCF_001546485_1_ASM154648v1_genomic:0.017705376,(GCF_000414445_1_ASM41444v1_genomic:0.000119495,GCF_000414465_1_ASM41446v1_genomic:0.000288031)1.000:0.018254326)1.000:0.187323158,GCF_003408845_1_ASM340884v1_genomic:0.273480559)1.000:0.031516103)1.000:0.105287618)1.000:0.184239398,(GCF_001546455_1_ASM154645v1_genomic:0.158834848,((((GCF_000414665_1_ASM41466v1_genomic:0.024400129,((GCF_002861905_1_ASM286190v1_genomic:0.019593979",
    ",GCF_013315115_1_ASM1331511v1_genomic:0.028631331)1.000:0.012662116,(GCF_000414585_1_ASM41458v1_genomic:0.000000005,GCF_000414605_1_ASM41460v1_genomic:0.000062682)1.000:0.024365016)0.997:0.009709274)1.000:0.017158240,(GCF_001546445_1_ASM154644v1_genomic:0.045083787,GCF_002861885_1_ASM286188v1_genomic:0.041601544)0.989:0.009811219)1.000:0.008444557,(GCF_000414625_1_ASM41462v1_genomic:0.042525597,GCF_003408835_1_ASM340883v1_genomic:0.049892270)0.999:0.010953734)1.000:0.036963531,((GCF_003426565_1_ASM342656v1_genomic:0.000997096,piotii_GCF_003397585_1_ASM339758v1_genomic:0.000227405)1.000:0.051433551,(GCF_000414545_1_ASM41454v1_genomic:0.049885899,((GCF_000414425_1_ASM41442v1_genomic:0.040267066,(GCF_000414485_1_ASM41448v1_genomic:0.041665521,GCF_000414505_1_ASM41450v1_genomic:0.029538413)1.000:0.011386125)0.275:0.004900178,(GCF_000414565_1_ASM41456v1_genomic:0.040629573,GCF_003397615_1_ASM339761v1_genomic:0.041114766)1.000:0.009480673)1.000:0.017839738)0.891:0.011092271)1.000:0.017665095)1.000:0.041690935)1.000:0.103183064)1.000:0.070463199",
    ",(GCF_000414705_1_ASM41470v1_genomic:0.049521408,(GCF_000414525_1_ASM41452v1_genomic:0.080724487,GCF_003408785_1_ASM340878v1_genomic:0.061182628)1.000:0.024820373)1.000:0.021054869)1.000:0.033807240)1.000:0.018980712,(GCF_000414645_1_ASM41464v1_genomic:0.027907653,GCF_013315075_1_ASM1331507v1_genomic:0.026610409)0.984:0.006296718)1.000:0.007287308,((((GCF_002862005_1_ASM286200v1_genomic:0.0,GCF_002862015_1_ASM286201v1_genomic:0.0):0.031702636,(GCF_003397605_1_ASM339760v1_genomic:0.018691990,GCF_013315045_1_ASM1331504v1_genomic:0.028742078)1.000:0.008607032)1.000:0.004053679,(GCF_003397665_1_ASM339766v1_genomic:0.023462907,(GCF_001660735_1_ASM166073v1_genomic:0.013671697,(GCF_003585755_1_ASM358575v1_genomic:0.022805192,(GCF_001660755_1_ASM166075v1_genomic:0.000031311,GCF_002861165_1_ASM286116v1_genomic:0.000000005)1.000:0.012753268)1.000:0.012810342)1.000:0.007322703)1.000:0.004506398)1.000:0.007041812,((GCF_002861965_1_ASM286196v1_genomic:0.0,GCF_002861975_1_ASM286197v1_genomic:0.0):0.027555977,GCF_003812765_1_ASM381276v1_genomic:0.015476481)1.000:0.005316359)0.983:0.002876082)0.313:0.002898759,((GCF_000414685_1_ASM41468v1_genomic:0.019432538,GCF_003408745_1_ASM340874v1_genomic:0.025614830)1.000:0.007383498,(GCF_000213955_1_ASM21395v1_genomic:0.021663753,GCF_013315025_1_ASM1331502v1_genomic:0.021564210)0.690:0.005300605)0.993:0.003197446)0.996:0.003120071);"
    )
  tree2 <- read.tree(text = text2)
  
  VisualizeMatching(JaccardRobinsonFoulds, tree1, tree2)

})

test_that("VisualizeMatching() works", {
  tree1 <- PectinateTree(1:11)
  tree2 <- tree1
  tree2$tip.label[c(11, 1)] <- tree1$tip.label[c(1, 11)]
  tree2r <- CollapseNode(tree2, 20:21)
  
  Minus <- function(...) {
    x <- MutualClusteringInfo(...)
    attr(x, "pairScores") <- -attr(x, "pairScores")
    x
  }
  expect_error(VisualizeMatching(Minus, PectinateTree(8), BalancedTree(8), 
                                 setPar = FALSE))
  
  skip_if_not_installed("vdiffr")
  skip_if(packageVersion("graphics") < "4.3")
  skip_if(packageVersion("vdiffr") < "1.0")
  
  TestVM <- function() {
    VisualizeMatching(MutualClusteringInfo, tree1, tree2, 
                      setPar = TRUE, precision = 3, matchZeros = FALSE,
                      Plot = plot.phylo)
  }
  vdiffr::expect_doppelganger("Test VM", TestVM)
  
  TestVMr <- function() {
    VisualizeMatching(MutualClusteringInfo, tree1, tree2r,
                      setPar = TRUE, precision = 3, matchZeros = TRUE, 
                      Plot = plot.phylo, cex = 1.5)
  }
  vdiffr::expect_doppelganger("Test VMr", TestVMr)
  
  vdiffr::expect_doppelganger("Visualize MCI matching", function() {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text="((1, 2), ((3, (4, 5)), (6, (7, (8, 9)))));")
    tree2 <- ape::read.tree(text="((1, 2), ((3, 4, (5, 9)), (6, (7, 8))));")
    VisualizeMatching(MutualClusteringInfo, tree1, tree2,
                      setPar = FALSE, precision = 3L,
                      Plot = TreeDistPlot,
                      matchZeros = FALSE,
                      leaveRoom = FALSE)
    VisualizeMatching(MutualClusteringInfo, tree2, tree1,
                      setPar = FALSE, precision = 3,
                      Plot = TreeDistPlot,
                      matchZeros = FALSE,
                      leaveRoom = FALSE)
  })
  
  vdiffr::expect_doppelganger("RF: Collapse a node", function() {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text="((1, 2), ((6, (7, 8)), (3, 4, (5, 9))));")
    tree2 <- ape::read.tree(text="((1, 2), ((3, (4, (5, 9))), (6, (7, 8))));")
    VisualizeMatching(RobinsonFouldsMatching, tree1, tree2,
                      setPar = FALSE, precision = 3,
                      Plot = TreeDistPlot,
                      matchZeros = TRUE,
                      leaveRoom = FALSE)
    VisualizeMatching(RobinsonFouldsMatching, tree2, tree1,
                      setPar = FALSE, precision = 3,
                      Plot = TreeDistPlot,
                      matchZeros = FALSE,
                      leaveRoom = FALSE)
  })
  
  
  vdiffr::expect_doppelganger("RF: Collapse and change", function() {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text="((1, 2), ((6, (7, 8)), (3, 4, (5, 9))));")
    tree2 <- ape::read.tree(text="((1, 2), ((3, (4, (5, 9))), ((6, 7), 8)));")
    VisualizeMatching(RobinsonFouldsMatching, tree1, tree2,
                      setPar = FALSE, precision = 3L,
                      Plot = TreeDistPlot,
                      matchZeros = TRUE,
                      leaveRoom = FALSE)
    VisualizeMatching(RobinsonFouldsMatching, tree2, tree1,
                      setPar = FALSE, precision = 3L,
                      Plot = TreeDistPlot,
                      matchZeros = FALSE,
                      leaveRoom = FALSE)
  })
  
  vdiffr::expect_doppelganger("RF VM Single splits; plainEdges", function() {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- ape::read.tree(text = "((1, 2), (3, 4, 5, 6, 7, 8));")
    tree2 <- ape::read.tree(text = "((1, 2, 3), (4, 5, 6, 7, 8));")
    VisualizeMatching(RobinsonFouldsMatching, tree1, tree2,
                      setPar = FALSE,
                      Plot = TreeDistPlot,
                      matchZeros = TRUE,
                      plainEdges = TRUE,
                      edge.width = NULL,
                      leaveRoom = FALSE)
    VisualizeMatching(RobinsonFouldsMatching, tree2, tree1,
                      setPar = FALSE,
                      Plot = TreeDistPlot,
                      matchZeros = FALSE,
                      plainEdges = FALSE,
                      leaveRoom = FALSE)
  })
  
  vdiffr::expect_doppelganger("JRF VM matchZeros FALSE", function() {
    JRF2 <- function(tree1, tree2, ...) 
      JaccardRobinsonFoulds(tree1, tree2, k = 2, allowConflict = FALSE, ...)
    
    tree1 <- RootTree(as.phylo(704564, 10), paste0("t", c(1, 4, 5, 8, 9)))
    tree2 <- RootTree(as.phylo(20165 , 10), paste0("t", c(1, 4)))
    VisualizeMatching(JRF2, tree1, tree2, matchZeros = FALSE)
  })
})

test_that("VisualizeMatching() handles unrooted trees", {
  skip_if_not_installed("graphics", "4.3")
  skip_if_not_installed("vdiffr", "1.0")
  
  vdiffr::expect_doppelganger("VM unrooted", function() {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- UnrootTree(BalancedTree(1:5))
    tree2 <- UnrootTree(PectinateTree(1:5))
    VisualizeMatching(RobinsonFouldsMatching, tree1, tree2,
                      setPar = FALSE,
                      Plot = TreeDistPlot)
  })
  
  vdiffr::expect_doppelganger("VM one rooted", function() {
    par(mfrow = c(2, 2), mar = rep(0.1, 4), cex = 1.5)
    tree1 <- UnrootTree(BalancedTree(1:5))
    tree2 <- PectinateTree(1:5)
    VisualizeMatching(RobinsonFouldsMatching, tree1, tree2,
                      setPar = FALSE,
                      Plot = TreeDistPlot)
  })
})
