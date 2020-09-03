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
load("SimulationDesign1.RData")

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

setwd("C:/Users/eddurso/Desktop/Projects/D'Urso MI and RS/Response Styles/Codes/Final Codes/DataSimulationStudy1")

for(i in 4:nrow(simdes)){
load(paste0("DataSim1_", i, ".RData"))
  
  #Run EFA up to a pre-determined number of factors
  EfaResults <- list()
  for(k in 1:length(Data$Data)){
    EfaResults[[k]] <- RunEFA(Data$Data[[k]], nfactors = 5)
  }
  
  if(simdes$scale[i] == "balanced"){
    balanced <-  T
  } else {
    balanced <- F
  }
  
  #Rotate results using various rotations (oblimin, target, partially specified target)
  RotatedResults <- list()
  RotatedResultsPoly <- list()
  for (k in 1:length(EfaResults)){
    RotatedResults[[k]] <-  rotating(EfaResults[[k]]$efaCor[[2]]$loadings, balanced = balanced)
    RotatedResultsPoly[[k]] <-  rotating(EfaResults[[k]]$efaPoly[[2]]$loadings, balanced = balanced)
  }
  
  save(EfaResults, RotatedResults, RotatedResultsPoly,
       file =paste0("C:/Users/eddurso/Desktop/Projects/D'Urso MI and RS/Response Styles/Codes/Final Codes/ResultsSimulationStudy1/",
                          "ResultsSim1_",i, ".RData"))  
  
    
}


