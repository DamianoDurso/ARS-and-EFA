##########################################################################################################3
#
#                     G E N E R A T I N G   D A T A    s I M S T U D Y  1  
#
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################source("SubRoutines.R")
setwd("C:/Users/eddurso/Desktop/Projects/D'Urso MI and RS/Response Styles/Codes/Final Codes")

library(mirt)
library(parallel)
source("SubRoutines.R")

#Factors simulation study
samplesize <- c(250, 500, 10000)
categories <- c(3,5,7)
scale <- c('balanced', 'unbalanced')
nitems <- c(12, 24)
categories <- c(3,5,7)
ARS <- c(0,0.38, 0.62, 1)

Simulation <- expand.grid(sample = samplesize, 
                          scale = scale, 
                          j = nitems, 
                          c = categories,
                          ARS = ARS
)
simdes <- cbind.data.frame(seed=seq(1:nrow(Simulation)),Simulation) #add the seed number




for(i in 1:nrow(simdes)){
  Data <-  GenSimulationStudy1(seed = simdes$seed[i], 
                                        sample = simdes$sample[i],
                                        scale = simdes$scale[i],
                                        nitems = simdes$j[i],
                                        cat = simdes$c[i],
                                        ARS = simdes$ARS[i],
                                        nrep = 10)

save(Data, file =paste0("C:/Users/eddurso/Desktop/Projects/D'Urso MI and RS/Response Styles/Codes/Final Codes/DataSimulationStudy1/",
  "DataSim1_",i, ".RData"))  
    
}



