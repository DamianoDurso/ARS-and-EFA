##########################################################################################################3
##########################################################################################################3
##########################################################################################################3
##########################################################################################################3
##########################################################################################################3
##########################################################################################################3
#                                F U N C T I O N S
##########################################################################################################3
##########################################################################################################3
##########################################################################################################3
##########################################################################################################3
##########################################################################################################3
##########################################################################################################3

myTryCatch <- function(expr) {
  warn <- err <- "none"                                                     # if no errors occur, let warnings and errors be "none"
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e                                                             # if an error occurs, save it in err
    }), warning=function(w) {
      warn <<- w                                                            # if a warning occurs, save it in warn
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)                                
}




##########################################################################################################3
#
#                     D A T A        G E N E R A T I O N
#
GenSimulationStudy2 <- function(seed, sample, scale, nitems, cat, ARS, nrep){
  
  #First Generate the matrix of factor loadings/discrimnation parameter  
  if(ARS == 0){
    if(nitems == 12){
      if(scale == "balanced"){
        a <- matrix(c(rep(c(1,0,-1,0),3), rep(c(0,1,0,-1),3)),
                    ncol=2)
      } else {
        a <- matrix(c(rep(c(1,0,1,0),3), rep(c(0,1,0,1),3)),  
                    ncol=2)
      }
      
    } else {
      if(scale == "balanced"){
        a <- matrix(c(rep(c(1,0,-1,0),6), rep(c(0,1,0,-1),6)), 
                    ncol=2)
      } else {
        a <- matrix(c(rep(c(1,0,1,0),6), rep(c(0,1,0,1),6)), 
                    ncol=2)
      }
    }
    
  }  else {
    if(nitems == 12){
      if(scale == "balanced"){
        a <- matrix(c(rep(c(1,0,-1,0),3), rep(c(0,1,0,-1),3), rep(ARS, nitems)),
                    ncol=3)
      } else {
        a <- matrix(c(rep(c(1,0,1,0),3), rep(c(0,1,0,1),3), rep(ARS, nitems)), 
                    ncol=3)
      }
      
    } else {
      if(scale == "balanced"){
        a <- matrix(c(rep(c(1,0,-1,0),6), rep(c(0,1,0,-1),6), rep(ARS, nitems)), 
                    ncol=3)
      } else {
        a <- matrix(c(rep(c(1,0,1,0),6), rep(c(0,1,0,1),6), rep(ARS, nitems)), 
                    ncol=3)
      }
    }  
}
  
  
  
  #Distance across categories 
  if(cat==3){
    distance <- 2
  } else if(cat==5){
    distance <- 1.25
  } else {
    distance <- 1.5
  }
  
  #Create thresholds
  diffs <- t(apply(matrix(runif(nitems*(cat-1), distance, distance),
                          nitems), 1, cumsum))
  diffs <- -(diffs - rowMeans(diffs)) #thresholds without shifts
  d <- diffs + c(seq(-1, 1, (2/11)))  #add shift to thresholds
  
  #  Distance between extreme categories
  d[1,1] - d[1,2]
  
  
  #create the factor scores
  if(ARS == 0){
    set.seed(seed)
    theta <- rmvnorm(sample, c(0,0), sigma = diag(2))  
  } else {
    set.seed(seed)
    thetaCT <- rmvnorm(sample, c(0,0), sigma = diag(2))  
    #ARS distribution censored
    thetaRS <- rnorm(sample, 0, 1)
    for(j in 1:length(thetaRS)){
      if (thetaRS[j] <= 0){
        thetaRS[j] <- 0
      } else {
        NULL
      }
    }
    #hist(thetaRS)
    theta <- cbind(thetaCT, thetaRS)
  }
  
  dat <- list()
  for (k in 1:nrep){
    #Create latent variables scores
    #normal distribution for the content trait
    #Generate data using MIRT
    dat[[k]] <- simdata(a, d, Theta = as.matrix(theta), itemtype = 'graded')
  }
  
  
  results <- list(dat, a, d)
  names(results) <- c('Data', "loadings", "tresholds")
  
  return(results)
}

##########################################################################################################3
#
#                     E S T I M A T I N G    E F A
#

#Run (categorical) exploratory factor analysis 
RunEFA <- function(X, nfactors = 5){
  efacor <- list()
  efapolychor <- list()

  #use also Parell Analysis
  efacor$pa <- fa.parallel(X, cor = "cor")
  efapolychor$pa <- fa.parallel.poly(X)
  
    
  for (i in 1:efacor$pa$nfact+1){
    efacor[[i]] <- fa(X,nfactors=i,n.obs=nrow(X), rotate = 'none', fm = "ml", cor = "cor") 
  }

  for (i in 1:efapolychor$pa$nfact+1){
  efapolychor[[i]] <- fa(X,nfactors=i,n.obs=nrow(X), rotate = 'none', fm = "ml", cor = "poly")
  }  
  
  results <- list(efacor, efapolychor)
  names(results) <- c("efaCor", "efaPoly")
  return(results)
}


##########################################################################################################3
#
#     R E F L E C T I N G    A N D   C O N V E R T I N G   R O T A T E D    E F A      S O L U T I O N S   
#


#Functions to reflect and convert the loadings with opposite sign
#function to refelect and convert loading matrix
reflectandconvmulti <- function(X, balanced ){
  nitem <- dim(X)[1]
  #Create reflection matrix
  if(dim(X)[2]<3){
    if(nitem == 12){
      if(balanced == T){
        keys <- matrix(c(rep(c(1,0,-1,0),3), rep(c(0,1,0,-1),3)), ncol=2)
      } else {
        keys <- matrix(c(rep(c(1,0,1,0),3), rep(c(0,1,0,1),3)), ncol=2)
      }
      
    } else {
      if(balanced == T){
        keys <- matrix(c(rep(c(1,0,-1,0),6), rep(c(0,1,0,-1),6)), ncol=2)
      } else {
        keys <- matrix(c(rep(c(1,0,1,0),6), rep(c(0,1,0,1),6)), ncol=2)
      }
    }
    
    
  } else {
    if(nitem == 12){
      if(balanced == T){
        keys <- matrix(c(rep(1, nitem),rep(c(1,0,-1,0),3), rep(c(0,1,0,-1),3)),
                       ncol=3)
      } else {
        keys <- matrix(c(rep(1, nitem), rep(c(1,0,1,0),3), rep(c(0,1,0,1),3)), 
                       ncol=3)
      }
      
    } else {
      if(balanced == T){
        keys <- matrix(c(rep(1, nitem),rep(c(1,0,-1,0),6), rep(c(0,1,0,-1),6)), 
                       ncol=3)
      } else {
        keys <- matrix(c(rep(1, nitem),rep(c(1,0,1,0),6), rep(c(0,1,0,1),6)), 
                       ncol=3)
      }
    }
    
  }
  

  for(i in 1:dim(X)[2]){
    for(k in 1:nitem){
      if(sign(X[k,i])==sign(keys[k,i])){
        NULL
      } else{
        X[k,i] <- -X[k,i]
      }
    }
  }
  return(X)
}



##########################################################################################################3
#
#                     R O T A T I N G     E F A      S O L U T I O N S   
#
#Function to rotate EFA solutions (oblimin, target, partially specified target using content trait(CT) and ARS)
rotatingnoARS <- function(X, balanced){
  nitems <- nrow(X)
  if(balanced == T){ #specify the target rotation differently based on whether the scale is balanced or not
    if(nitems == 12){
      keystarg <- matrix(c(rep(c(1,0,-1,0),3), rep(c(0,1,0,-1),3)), ncol=2)
      keystargCT <- matrix(c(rep(c(NA,0,NA,0),3), rep(c(0,NA,0,NA),3)), ncol=2)

      obli <- oblimin(X)  
      Target <- myTryCatch(target.rot(X, keystarg))
      TargCT <- myTryCatch(targetQ(X, Target = keystargCT))

    } else {
      keystarg <- matrix(c(rep(c(1,0,-1,0),6), rep(c(0,1,0,-1),6)), ncol=2)
      keystargCT <- matrix(c(rep(c(NA,0,NA,0),6), rep(c(0,NA,0,NA),6)), ncol=2)

      obli <- oblimin(X)  
      Target <- myTryCatch(target.rot(X, keystarg))
      TargCT <- myTryCatch(targetQ(X, Target = keystargCT))

    }
  } else {
    if(nitems == 12){
      keystarg <- matrix(c(rep(c(1,0,1,0),3), rep(c(0,1,0,1),3)), ncol=2)
      keystargCT <- matrix(c(rep(c(NA,0,NA,0),3), rep(c(0,NA,0,NA),3)), ncol=2)

      obli <- oblimin(X)  
      Target <- myTryCatch(target.rot(X, keystarg))
      TargCT <- myTryCatch(targetQ(X, Target = keystargCT))

    } else {
      keystarg <- matrix(c(rep(c(1,0,1,0),6), rep(c(0,1,0,1),6)), ncol=2)
      keystargCT <- matrix(c(rep(c(NA,0,NA,0),6), rep(c(0,NA,0,NA),6)), ncol=2)

      obli <- oblimin(X)  
      Target <- myTryCatch(target.rot(X, keystarg))
      TargCT <- myTryCatch(targetQ(X, Target = keystargCT))

    }
    
  }
  
  
  #Functions to reflect and convert the loadings with opposite sign
#  oblirefl <- reflectandconvmulti(obli$loadings, balanced)
  
  if(Target$error == "none"){
    Targrefl <- reflectandconvmulti(Target$value$loadings, balanced)
  } else {
    Targrefl <- "error found" 
  }
  
  if(TargCT$error == "none"){
    TargCTrefl <- reflectandconvmulti(TargCT$value$loadings, balanced)
  } else {
    TargCTrefl <- "error found" 
  }
  
  
  results <- list(obli, Target, TargCT, Targrefl, TargCTrefl)
  names(results) <- c("Oblimin", "Target", "TargetCT", 
                      "TargetRef", "TargetCTRef")
  return(results)
  
} 


#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#


#Function to rotate EFA solutions (oblimin, target, partially specified target using content trait(CT) and ARS)
rotatingARS <- function(X, balanced){
  nitems <- nrow(X)
  if(balanced == T){ #specify the target rotation differently based on whether the scale is balanced or not
    if(nitems == 12){
      keystarg <- matrix(c(rep(1, nitems), rep(c(1,0,-1,0),3), rep(c(0,1,0,-1),3)), ncol=3)
      keystargCT <- matrix(c(rep(NA, nitems),rep(c(NA,0,NA,0),3), rep(c(0,NA,0,NA),3)), ncol=3)
      keystargCTARS <- matrix(c(rep(1, nitems),rep(c(NA,0,NA,0),3), rep(c(0,NA,0,NA),3)), ncol=3)
      keystargARS <- matrix(c(rep(1, nitems),rep(c(NA,NA,NA,NA),3), rep(c(NA,NA,NA,NA),3)), ncol=3)
      
      obli <- oblimin(X)  
      Target <- myTryCatch(target.rot(X, keystarg))
      TargCT <- myTryCatch(targetQ(X, Target = keystargCT))
      TargCTARS <- myTryCatch(targetQ(X, Target = keystargCTARS))
      TargARS <- myTryCatch(targetQ(X, Target = keystargARS))
      
      } else {
        keystarg <- matrix(c(rep(1, nitems), rep(c(1,0,-1,0),6), rep(c(0,1,0,-1),6)), ncol=3)
        keystargCT <- matrix(c(rep(NA, nitems),rep(c(NA,0,NA,0),6), rep(c(0,NA,0,NA),6)), ncol=3)
        keystargCTARS <- matrix(c(rep(1, nitems),rep(c(NA,0,NA,0),6), rep(c(0,NA,0,NA),6)), ncol=3)
        keystargARS <- matrix(c(rep(1, nitems),rep(c(NA,NA,NA,NA),6), rep(c(NA,NA,NA,NA),6)), ncol=3)
        
        obli <- oblimin(X)  
        Target <- myTryCatch(target.rot(X, keystarg))
        TargCT <- myTryCatch(targetQ(X, Target = keystargCT))
        TargCTARS <- myTryCatch(targetQ(X, Target = keystargCTARS))
        TargARS <- myTryCatch(targetQ(X, Target = keystargARS))
        
      }
  } else {
    if(nitems == 12){
      keystarg <- matrix(c(rep(1, nitems), rep(c(1,0,1,0),3), rep(c(0,1,0,1),3)), ncol=3)
      keystargCT <- matrix(c(rep(NA, nitems),rep(c(NA,0,NA,0),3), rep(c(0,NA,0,NA),3)), ncol=3)
      keystargCTARS <- matrix(c(rep(1, nitems),rep(c(NA,0,NA,0),3), rep(c(0,NA,0,NA),3)), ncol=3)
      keystargARS <- matrix(c(rep(1, nitems),rep(c(NA,NA,NA,NA),3), rep(c(NA,NA,NA,NA),3)), ncol=3)
      
      obli <- oblimin(X)  
      Target <- myTryCatch(target.rot(X, keystarg))
      TargCT <- myTryCatch(targetQ(X, Target = keystargCT))
      TargCTARS <- myTryCatch(targetQ(X, Target = keystargCTARS))
      TargARS <- myTryCatch(targetQ(X, Target = keystargARS))
      
    } else {
      keystarg <- matrix(c(rep(1, nitems), rep(c(1,0,1,0),6), rep(c(0,1,0,1),6)), ncol=3)
      keystargCT <- matrix(c(rep(NA, nitems),rep(c(NA,0,NA,0),6), rep(c(0,NA,0,NA),6)), ncol=3)
      keystargCTARS <- matrix(c(rep(1, nitems),rep(c(NA,0,NA,0),6), rep(c(0,NA,0,NA),6)), ncol=3)
      keystargARS <- matrix(c(rep(1, nitems),rep(c(NA,NA,NA,NA),6), rep(c(NA,NA,NA,NA),6)), ncol=3)
      
      obli <- oblimin(X)  
      Target <- myTryCatch(target.rot(X, keystarg))
      TargCT <- myTryCatch(targetQ(X, Target = keystargCT))
      TargCTARS <- myTryCatch(targetQ(X, Target = keystargCTARS))
      TargARS <- myTryCatch(targetQ(X, Target = keystargARS))
      
    }
  
  }
    

    #Functions to reflect and convert the loadings with opposite sign
#    oblirefl <- reflectandconvmulti(obli$loadings, balanced)
    
    if(Target$error == "none"){
      Targrefl <- reflectandconvmulti(Target$value$loadings, balanced = balanced)
    } else {
      Targrefl <- "error found" 
    }
    
    if(TargCT$error == "none"){
      TargCTrefl <- reflectandconvmulti(TargCT$value$loadings, balanced = balanced)
    } else {
      TargCTrefl <- "error found" 
    }

   if(TargCTARS$error == "none"){
   TargCTARSrefl <- reflectandconvmulti(TargCTARS$value$loadings, balanced = balanced)
   } else {
     TargCTARSrefl <- "error found" 
   }
  
   if(TargARS$error == "none"){
      TargARSrefl <- reflectandconvmulti(TargARS$value$loadings, balanced = balanced)
    } else {
      TargARSrefl <- "error found" 
    }
    
    
    results <- list(obli, Target, TargCT, TargCTARS, TargARS, Targrefl, TargCTrefl,TargCTARSrefl, TargARSrefl)
    names(results) <- c("Oblimin", "Target", "TargetCT","TargetCTARS", "TargetARS","TargetRef", 
                        "TargetCTRef","TargetCTARSRef", "TargetARSRef")
    return(results)
    
  } 
##########################################################################################################3
##########################################################################################################3


#USE THE CHULL METHOD WITH KMO AS PROPOSED BY LORENZO SEVA AND TIMMERMAN (2012)
#function to calculate chull
ChullCAF <- function(X, maxfac = 5){
  
  var <- ncol(X[[1]]$r)
  factors <- seq(1:maxfac)
  
  Chullmat <- matrix(NA, ncol = 2, nrow = maxfac+1)
  colnames(Chullmat) <- c("dof", "CAF")
  
  #calculate DOF for each model (formula from Lorenz-seva 2011)
  for(i in 0:maxfac){
    Chullmat[i+1,1] <- var*i - 1/2*i*(i-1)
  }
  
  #calculate CAF [COMMON PART ACCOUNTED FOR] for NULL model
  Chullmat[1,2] <- 1-KMO(X[[1]]$r)$MSA
  
  #create a list with efa solution
  efamods <- list()
  for(i in 1:maxfac){
    efamods[[i]] <- X[[i]] 
  }
  
  #calculate CHull
  for(i in 1:maxfac+1){
    Chullmat[i,2] <- 1-KMO(efamods[[i-1]]$residual)$MSA
  }
  #___________________________________________________________________________________________________________#
  #___________________________________________________________________________________________________________#
  #___________________________________________________________________________________________________________#
  
  #  ###Now the same but for ordinal data
  #  CHullmatcat <- matrix(NA, ncol = 2, nrow = maxfac+1)
  #  CHullmatcat[,1] <- Chullmat[,1]
  #  colnames(CHullmatcat) <- c("dof", "CAF")
  #  CHullmatcat[1,2] <- 1-KMO(X)$MSA
  
  #  efamodscat <- list()
  #  for(i in 1:maxfac){
  #    efamodscat[[i]] <- X[[i]] 
  #  }
  
  #  #calculate CHull
  #  for(i in 2:5){
  #    CHullmatcat[i,2] <- 1-KMO(efamodscat[[i-1]]$residual)$MSA
  #  }
  
  
  #chull obtained
  outchull <- CHull(Chullmat, bound = "upper")
  #  plot(outchull)
  
  #  outchullcat <- CHull(CHullmatcat, bound = "upper")
  #  plot(outchull)
  
  
  #indicate selected solution (-1 because 1 is null model)
  if(is.list(outchull)){
    solution <- which(Chullmat[,1] == max(outchull$Solution$complexity))-1
    #  solutioncat <- which(CHullmatcat[,1] == max(outchullcat$Solution$complexity))-1
    out <- c(outchull, solution)
    names(out)[6] <- c("selectedmod")
  } else {
    out <- list()
    out[[1]] <- outchull
    out[[2]] <- "none"
    names(out) <- c("error", "selectedmod")
  }
  return(out)
}
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

