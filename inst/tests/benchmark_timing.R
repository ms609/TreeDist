if (FALSE){  ## Don't allow codecov to get drawn in to this file
  load_all()
data('inapplicable.datasets')
scores <- c(
"Agnarsson2004" =  778 , # 0.044 mins   
"Aria2015" =       143 , # 0.041        
"Asher2005" =      345 , # 0.00328      
"Capa2011" =       385 , # 0.131        
"Conrad2008" =     1761, # 0.119        
"DeAssis2011" =    64  , # 0            
"Aguado2009" =     579 , # 4.979        
"Dikow2009" =      1611, # 22.35        # This is a slow one [22 mins]
"Eklund2004" =     440 , # 0.080        
"Geisler2001" =    1295, # 115.3        # This is a slow one [115 mins]
"Giles2015" =      710 , # 9.80         
"Griswold1999" =   407 , # 0.020        
"Liljeblad2008" =  2868, # 1.768        
"Loconte1991" =    539 , # 0.208        
"Longrich2010" =   131 , # 0.0033       
"OLeary1999" =     508 , # 0.158        
"OMeara2014" =     273 , #              
"Rougier2012" =    1215, #              
"Rousset2004" =    259 , #              
"Sano2011" =       223 , #              
"Sansom2010" =     189 , #              
"Schulze2007" =    164 , #              
"Shultz2007" =     454 , #              
"Vinther2008" =    79  , #              
"Wetterer2000" =   559 , #              
"Wills2012" =      273 , #              
"Wilson2003" =     879 , #              
"Wortley2006" =    482 , #              
"Zanol2014" =      1311, #              
"Zhu2013" =        638 ) #              

nj.tree <- lapply(inapplicable.phyData, NJTree)
timestart <- double(length(scores))
timeend   <- double(length(scores))
names(timestart) <- names(timeend) <- names(scores)

#install_github('ms609/inapplicable', rel='cefb5669352aca6425516805f60108063383b6c2')

## Find a good searchHits value
dataName <- names(scores)[1]

candidates <- c(4, 6, 8, 10, 12, 15, 18, 22, 26, 30, 35, 40, 50, 60, 75, 90, 120, 150)
bench <- function (searchHits, bootHits = searchHits) system.time(Ratchet(nj.tree[[dataName]], inapplicable.phyData[[dataName]],                    swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap), stopAtScore=scores[[dataName]], ratchHits=1000, ratchIter=10000, searchIter=3200, 
                                                              searchHits=searchHits, bootstrapHits=bootHits, verbosity=1L))
bench2 <- function (bootHits) system.time(Ratchet(nj.tree[[dataName]], inapplicable.phyData[[dataName]],                    swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap), stopAtScore=scores[[dataName]], ratchHits=1000, ratchIter=10000, searchIter=3200, 
                                                  searchHits=30, bootstrapHits=bootHits, verbosity=1L))


results1 <- vapply(candidates, bench, double(5))[1, ]
results2 <- vapply(candidates, bench2, double(5))[1, ]
results3.10 <- vapply(candidates, bench, bootHits = 10, double(5))[1, ]
results3.20 <- vapply(candidates, bench, bootHits = 20, double(5))[1, ]
results3.30 <- vapply(candidates, bench, bootHits = 30, double(5))[1, ]
results3.50 <- vapply(candidates, bench, bootHits = 50, double(5))[1, ]
results3.90 <- vapply(candidates, bench, bootHits = 90, double(5))[1, ]

greens <- c('#00ddbb', '#33dd99', '#99dd33', '#00ddbb', '#22bb22')
plot(results1 ~ candidates, ylab='Time taken / s', col=2, pch=4)
points(results2 ~ candidates, pch=3, col=4)
text(labels="10", results3.10 ~ candidates, pch=6, col=greens[1])
text(labels="20", results3.20 ~ candidates, pch=6, col=greens[2])
text(labels="30", results3.30 ~ candidates, pch=6, col=greens[3])
text(labels="50", results3.50 ~ candidates, pch=6, col=greens[4])
text(labels="90", results3.90 ~ candidates, pch=6, col=greens[5])
legend('topright', legend=c('SearchHits (bh=sh)', 'BootHits (sh=30)',
                            paste0('SearchHits (bh=', c(10, 20, 30, 50, 90), ')')), 
       pch=c(4, 3, rep(6, 5)), 
       col=c(2,4, greens))

n <- length(candidates)
allRes <- array(c(results1, results2, results3.10, results3.20, results3.30, results3.50, results3.90,
               candidates, rep(30, n),  candidates, candidates, candidates, candidates, candidates, 
               candidates, candidates, rep(10, n), rep(20, n), rep(30, n), rep(50, n), rep(90, n)),
             c(7*length(candidates), 3)) 

arr <- allRes
rgl::plot3d(arr[, 2], arr[, 3], arr[, 1], zlab='Time', xlab='SearchHits', ylab='BootHits')
manyBootHits <- allRes[allRes[, 3] > 30, ]
arr <- manyBootHits
rgl::plot3d(arr[, 2], arr[, 3], arr[, 1], zlab='Time', xlab='SearchHits', ylab='BootHits')

library('profvis')
RRprofStart()
Rprof()
oTree <- Ratchet(nj.tree[[dataName]], inapplicable.phyData[[dataName]],
                 swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
                 stopAtScore=scores[[dataName]],
                 ratchHits=1000, ratchIter=10000,
                 searchIter=3200, searchHits=12, verbosity=2L)
Rprof(NULL)
summaryRprof()
RRprofStop()
RRprofReport()

pd <- proftools::profileExpr(
oTree <- Ratchet(nj.tree[[dataName]], inapplicable.phyData[[dataName]],
                 swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
                 stopAtScore=scores[[dataName]],
                 ratchHits=1000, ratchIter=10000,
                 searchIter=3200, searchHits=12, verbosity=2L)
)
dev.new()
library(proftools)
proftools::flameGraph(pd)
proftools::plotProfileCallGraph(pd, style=google.style, score='total')

proftools::plotProfileCallGraph(pd)
proftable

for (dataName in names(scores)) {
 cat("\n ========", format(Sys.time(), "%b %d %X"), ":", dataName, ": Target", scores[dataName], "========\n")
 timestart[dataName] <- Sys.time()
 oTree <- Ratchet(nj.tree[[dataName]], inapplicable.phyData[[dataName]],
                  swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
                  stopAtScore=scores[[dataName]],
                  ratchHits=1000, ratchIter=10000,
                  searchIter=3200, searchHits=12, verbosity=2L)
 timeend[dataName] <- Sys.time()
 cat("\n > Time taken: ", (timeend[dataName] - timestart[dataName]) / 60, "mins\n")
}

## TESTING
{
dataset <- names(scores)[4]
 cat("\n ========", format(Sys.time(), "%b %d %X"), ":", dataset, "========\n")
 timestart[dataset] <- Sys.time()
 oTree <- Ratchet(nj.tree[[dataset]], inapplicable.phyData[[dataset]],
                  swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
                  stopAtScore=scores[[dataset]], 
                  ratchHits=1000, ratchIter=10000,
                  searchIter=3200, searchHits=12, verbosity=4L)
 timeend[dataset] <- Sys.time()
 cat("\n > Time taken: ", (timeend[dataset] - timestart[dataset]) / 60, "mins\n")
}

} # End top-level IF