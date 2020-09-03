##########################################################################################################3
#
#                     G E N E R A T I N G   D A T A    s I M S T U D Y  2 
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
library(MASS)
library(mvtnorm)
source("SubRoutines2.R")


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
simdesmulti <- cbind.data.frame(seed=seq(1:nrow(Simulation)),Simulation) #add the seed number
#save(simdesmulti, file = "SimulationDesign2.RData")



for(i in 1:nrow(simdesmulti)){
  Data <-  GenSimulationStudy2(seed = simdesmulti$seed[i], 
                               sample = simdesmulti$sample[i],
                               scale = simdesmulti$scale[i],
                               nitems = simdesmulti$j[i],
                               cat = simdesmulti$c[i],
                               ARS = simdesmulti$ARS[i],
                               nrep = 10)
  
  save(Data, file =paste0("C:/Users/eddurso/Desktop/Projects/D'Urso MI and RS/Response Styles/Codes/Final Codes/DataSimulationStudy2/",
                          "DataSim2_",i, ".RData"))  
  
}



