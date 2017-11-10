data(inapplicable.datasets)
scores <- c(
"Agnarsson2004" =  778 , # 0.044 mins             
"Aguado2009" =     579 , # 4.979             
"Aria2015" =       143 , # 0.041             
"Asher2005" =      345 , # 0.00328      
"Capa2011" =       385 , # 0.131       
"Conrad2008" =     1761, # 0.119        
"DeAssis2011" =    64  , # 0            
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

nj.tree <- lapply(inapplicable.phyData, TreeSearch::NJTree)
timestart <- double(length(scores))
timeend   <- double(length(scores))
names(timestart) <- names(timeend) <- names(scores)

install_github('ms609/inapplicable', rel='cefb5669352aca6425516805f60108063383b6c2')

for (dataset in names(inapplicable.phyData)) {
  cat("\n ========", format(Sys.time(), "%b %d %X"), ":", dataset, "========\n")
  timestart[dataset] <- Sys.time()
  oTree <- RatchetSearch(nj.tree[[dataset]], inapplicable.phyData[[dataset]], stopAtScore=scores[[dataset]],
  k=1000, maxIt=10000, maxIter=3200, maxHits=12, verbosity=0)
  timeend[dataset] <- Sys.time()
  cat("\n > Time taken: ", (timeend[dataset] - timestart[dataset]) / 60, "mins\n")
}

timetaken <- timeend - timestart
sum(timetaken)

                   