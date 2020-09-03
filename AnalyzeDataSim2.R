##########################################################################################################3
#
#                     E S T I M A T I N G    E F A   A N D   R O T A T I N G   
#
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
setwd("C:/Users/eddurso/Desktop/Projects/D'Urso MI and RS/Response Styles/Codes/Final Codes")
source("SubRoutines.R")
library(parallel)
library(GPArotation)
library(psych)
load("SimulationDesign2.RData")

#Factors simulation study
#samplesize <- c(250, 500, 10000)
#categories <- c(3,5,7)
#scale <- c('balanced', 'unbalanced')
#nitems <- c(12, 24)
#categories <- c(3,5,7)
#ARS <- c(0,0.38, 0.62, 1)

#Simulation <- expand.grid(sample = samplesize, 
#                          scale = scale, 
#                          j = nitems, 
#                          c = categories,
#                          ARS = ARS
#)
#simdes <- cbind.data.frame(seed=seq(1:nrow(Simulation)),Simulation) #add the seed number
#save(simdes, file = "SimulationDesign1.RData")

setwd("C:/Users/eddurso/Desktop/Projects/D'Urso MI and RS/Response Styles/Codes/Final Codes/DataSimulationStudy2")

for(i in 4:nrow(simdesmulti)){
  load(paste0("DataSim2_", i, ".RData"))
  
  #Run EFA up to a pre-determined number of factors
  tic()
  EfaResults <- list()
  for(k in 1:length(Data$Data)){
    EfaResults[[k]] <- RunEFA(Data$Data[[k]], nfactors = 5)
  }
  time <- toc()
  
  if(simdesmulti$scale[i] == "balanced"){
    balanced <-  T
  } else {
    balanced <- F
  }
  
  #Rotate results using various rotations (oblimin, target, partially specified target)
  RotatedResults <- list()
  RotatedResultsPoly <- list()
  NOARSRotatedResults <- list()
  NOARSRotatedResultsPoly <- list()
  
  for (k in 1:length(EfaResults)){
    RotatedResults[[k]] <-  rotatingARS(EfaResults[[k]]$efaCor[[3]]$loadings, balanced = balanced)
    RotatedResultsPoly[[k]] <-  rotatingARS(EfaResults[[k]]$efaPoly[[3]]$loadings, balanced = balanced)
    NOARSRotatedResults[[k]] <- rotatingnoARS(EfaResults[[k]]$efaCor[[2]]$loadings, balanced = balanced)
    NOARSRotatedResultsPoly[[k]] <- rotatingnoARS(EfaResults[[k]]$efaPoly[[2]]$loadings, balanced = balanced)
  }

  
  save(EfaResults, RotatedResults, RotatedResultsPoly, NOARSRotatedResults, NOARSRotatedResultsPoly,
       file =paste0("C:/Users/eddurso/Desktop/Projects/D'Urso MI and RS/Response Styles/Codes/Final Codes/ResultsSimulationStudy2/",
                    "ResultsSim2_",i, ".RData"))  
  
  
}


