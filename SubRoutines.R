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
GenSimulationStudy1 <- function(seed, sample, scale, nitems, cat, ARS, nrep){
  
  #First Generate the matrix of factor loadings/discrimnation parameter  
  if(ARS == 0){
    if(scale =="balanced"){
      a <- matrix(c(rep(1,nitems/2), rep(-1,nitems/2)), ncol = 1, nrow = nitems)
    } else {
      a <- matrix(c(rep(1,nitems)), ncol = 1, nrow = nitems)
    }
    
  }  else {
    if(scale == "balanced"){
      a <- matrix(c(rep(1,nitems/2), rep(-1, nitems/2), 
                    rep(ARS,nitems)), 
                  ncol = 2, nrow = nitems)
    } else {
      a <- matrix(c(rep(1,nitems), 
                    rep(ARS, nitems)), 
                  ncol = 2, nrow = nitems)
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
    theta <- rnorm(sample, 0, 1)  
  } else {
    set.seed(seed)
    thetaCT <- rnorm(sample, 0, 1)  
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

remove <- NULL    
 for(b in 1:length(dat)){
   if(all(apply(dat[[b]], 2, function(a) length(unique(a))==length(unique(dat[[b]][,1]))))){
     NULL
   } else {
     remove[b] <- b 
   }
 }

 if(is.null(remove)){
   NULL
 } else {
   nas <- which(is.na(remove))
   remove <- remove[-c(nas)]
   dat <- dat[-c(remove)]
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
  
  for (i in 1:nfactors){
    efacor[[i]] <- fa(X,nfactors=i,n.obs=nrow(X), rotate = 'none', fm = "ml", cor = "cor") 
    efapolychor[[i]] <- fa(X,nfactors=i,n.obs=nrow(X), rotate = 'none', fm = "ml", cor = "poly")
  }
  
  #use also Parell Analysis
  efacor$pa <- fa.parallel(X, cor = "cor")
  efapolychor$pa <- fa.parallel(X, cor = "poly")
  
  
  results <- list(efacor, efapolychor)
  names(results) <- c("efaCor", "efaPoly")
  return(results)
}


##########################################################################################################3
#
#     R E F L E C T I N G    A N D   C O N V E R T I N G   R O T A T E D    E F A      S O L U T I O N S   
#


#Functions to reflect and convert the loadings with opposite sign
reflectandconv <- function(X){
  nitem <- dim(X)[1]
  #Create reflection matrix
  if(dim(X)[2]>1){
    keys <- matrix(c(rep(1,nitem/2), rep(-1,nitem/2), rep(1,nitem)), ncol=2)
  } else {
    keys <- matrix(c(rep(1,nitem/2), rep(-1, nitem/2), ncol=1))
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
rotating <- function(X, balanced){
  items <- nrow(X)
  if(balanced == T){ #specify the target rotation differently based on whether the scale is balanced or not
    keystarg <- matrix(c(rep(1, items/2), rep(-1, items/2), rep(1, items)), ncol= 2, nrow = items)
    keystargCT <-  matrix(c(rep(1, items/2), rep(-1, items/2), rep(NA, items)), ncol= 2, nrow = items)
    keystargARS <- matrix(c(rep(NA, items), rep(1, items)), ncol= 2, nrow = items)
    
    obli <- oblimin(X)  
    Target <- myTryCatch(target.rot(X, keystarg))
    TargCT <- myTryCatch(targetQ(X, Target = keystargCT))
    TargARS <- myTryCatch(targetQ(X, Target = keystargARS))
    
    #Functions to reflect and convert the loadings with opposite sign
    
    
    oblirefl <- reflectandconv(obli$loadings)
    
    if(Target$error == "none"){
      Targrefl <- reflectandconv(Target$value$loadings)
    } else {
      Targrefl <- "error found" 
    }
    
    if(TargCT$error == "none"){
      TargCTrefl <- reflectandconv(TargCT$value$loadings)
    } else {
      TargCTrefl <- "error found" 
    }
    
    if(TargARS$error == "none"){
      TargARSrefl <- reflectandconv(TargARS$value$loadings)
    } else {
      TargARSrefl <- "error found" 
    }
    
    
    results <- list(obli, Target, TargCT, TargARS, oblirefl, Targrefl, TargCTrefl, TargARSrefl)
    names(results) <- c("Oblimin", "Target", "TargetCT", "TargetARS", 
                        "ObliminRef", "TargetRef", "TargetCTRef", "TargetARSRef")
    return(results)
    
  } else {
    
#    keystarg <- matrix(c(rep(1, items), rep(1, items)), ncol= 2, nrow = items)
    keystargCT <-  matrix(c(rep(1, items), rep(NA, items)), ncol= 2, nrow = items)
    keystargARS <- matrix(c(rep(NA, items), rep(1, items)), ncol= 2, nrow = items)
    
    obli <- oblimin(X)  
#    Target <- myTryCatch(target.rot(X, keystarg))
    TargCT <- myTryCatch(targetQ(X, Target = keystargCT))
    TargARS <- myTryCatch(targetQ(X, Target = keystargARS))
    
    #Functions to reflect and convert the loadings with opposite sign
    
    
    oblirefl <- reflectandconv(obli$loadings)
    Targetrefl <- reflectandconv(X)
    
#    if(Target$error == "none"){
#      Targrefl <- reflectandconv(Target$value$loadings)
#    } else {
#      Targrefl <- "error found" 
#    }
    
    if(TargCT$error == "none"){
      TargCTrefl <- reflectandconv(TargCT$value$loadings)
    } else {
      TargCTrefl <- "error found" 
    }
    
    if(TargARS$error == "none"){
      TargARSrefl <- reflectandconv(TargARS$value$loadings)
    } else {
      TargARSrefl <- "error found" 
    }
    
    
    results <- list(obli, TargCT, TargARS, oblirefl, TargCTrefl, TargARSrefl)
    names(results) <- c("Oblimin", "TargetCT", "TargetARS", 
                        "ObliminRef", "TargetCTRef", "TargetARSRef")
    return(results)
    
  }
  
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

##########################################################################################################3
##########################################################################################################3
#Function to rotate and converti uni and multidimensional loadings without ARS

reflectandconv <- function(X){
  nitem <- dim(X)[1]
  #Create reflection matrix
  if(dim(X)[2]>1){
    keys <- matrix(c(rep(1,nitem/2), rep(-1,nitem/2), rep(1,nitem)), ncol=2)
  } else {
    keys <- matrix(c(rep(1,nitem/2), rep(-1, nitem/2), ncol=1))
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



########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

