# scenario set up
library(zoo)
library(abind)
library(dplyr)
#------------------------------------------------------------------------------#
# Functions for building prey/wound scenarios----
#------------------------------------------------------------------------------#

# function for building prey scenarios
build.scenario.prey <- function(
  preyDat, #plankton indices
  ref_yrs = list(2010:2019), #reference years
  proj_yrs = list(1:100), # the projected years for which ref_yrs apply
  preyChange = 1 # relative change in prey availability
  ){
  print(paste0("Number of years projected, nT = ",max(unlist(proj_yrs))))
  #browser()
  if (any(!sapply(c("nT","nBoot"),exists))){
    stop("!!! projection specs missing?")}
  if (any(!1:nT %in% unlist(proj_yrs))){
    warning("CHECK: some projection years (1:nT) unspecified?")}
  if (length(ref_yrs) != length(proj_yrs)){
    stop("!!! projection year periods do not correspond to reference year periods")}
  if (any(unlist(ref_yrs) <1980)){
    stop("!!! ref_yrs need to be 4-digit years (e.g., 1999)")}
  if (!exists("seed")){
    seed <- 2022
    warning(paste0("seed missing, using seed=",seed))
  }
  
  # mean/sd for standardizing (used in model estimation)
  preyDat_roll <- data.frame(
    year = preyDat$year,
    missing = preyDat$missing,
    preyDat %>% select(-year,-missing) %>%
      rollapply(3,function(x){mean(x,na.rm=T)},partial=T,align="right")) %>% 
    filter(year %in% 1990:2019)
  foodMean <- apply(preyDat_roll[,c("food1","food2")],2,mean)
  foodStdv <- apply(preyDat_roll[,c("food1","food2")],2,sd)

  # raw prey values sampled
  preyArrayList <- list()
  
  set.seed(seed) #ensure we can replicate
  for (i in 1:length(ref_yrs)){
    if (sum(preyDat$missing[which(preyDat$year %in% ref_yrs[[i]])]) == length(ref_yrs[[i]])){
      stop(paste0("!!! all ref_yrs in period ",i, " have missing prey index values"))}
    Yrs <- proj_yrs[[i]]; nYrs <- length(Yrs)
    if (i==1){
      preyArrayList[[i]] <- array(
        # values (need to add 2 at year=1 for the 3-yr rolling average)
        # (rep is needed to make sample fn work if ref_yrs is a single year)
        preyDat[sample(rep(which(preyDat$year %in% ref_yrs[[i]] & preyDat$missing==0
                             ),2), nBoot * (nYrs + 2), replace = TRUE), c("food1","food2")] %>% as.matrix(),
        # dimensions
        c(nYrs + 2, nBoot, 2),
        # names
        list(c(NA, NA, Yrs), NULL, c("food1", "food2")))
    } else {
      preyArrayList[[i]] <- array(
        # values
        # (rep is needed to make sample fn work if ref_yrs is a single year)
        preyDat[sample(rep(which(preyDat$year %in% ref_yrs[[i]] & preyDat$missing==0
                             ),2), nBoot * (nYrs), replace = TRUE), c("food1","food2")] %>% as.matrix(),
        # dimensions
        c(nYrs, nBoot, 2),
        # names
        list(c(Yrs), NULL, c("food1", "food2")))
    }
  }
  # create the raw index array, and change by some quantity (e.g., noise effects)
  preyArray <- abind(preyArrayList,along=1) + log(preyChange)
  
  # calculate the rolling 3-yr average (previous 3 years, hence "right")
  preyArray_roll <- apply(preyArray, c(2,3), rollmean, k = 3, align = "right")

  # standardize by the values used in the model estimation
  preyArray_roll <- sweep(preyArray_roll, 3, c(foodMean[1], foodMean[2]), "-")
  preyArray_roll <- sweep(preyArray_roll, 3, c(foodStdv[1], foodStdv[2]), "/")
  preyArray_perm <- aperm(preyArray_roll, c(2,1,3))
  
  return(round(preyArray_perm,4))
}


# function for building wound scenarios
build.scenario.wound <- function(
  iTheta,  #posterior mortality rates from model estimation
  eps.i,
  proj_yrs = list(1:100), # the projected years for which change lists apply
  iE.change = list(c(1,1,1,1,1)), # change in state-specific entanglement injury rate
  iV.change = list(c(1,1,1,1,1)), # change in state-specific vessel-strike injury rate
  iV.change.t = list(c(1,1,1,1,1)) # annual change in state-specific vessel-strike injury rate
){
  #assumes some projection specs have already been defined
  if (any(!sapply(c("nT","nBoot","nWoundStates", "nEntangleStages"),exists))){
    stop("!!! some projection specs missing?")}
  # if (any(!sapply(c("eps.i"),exists))){
  #   stop("!!! random temporal variation for injury rates missing")}
  
  woundArrayList <- list()
  
  for (period in 1:length(proj_yrs)){
    woundArrayList[[period]] <- array(NA, c(nBoot, nWoundStates, nWoundStates, nEntangleStages, nT), 
                                      list(NULL, woundStates, woundStates, entangleStages, year=1:nT))
    
    for (i in 1:nBoot){
      for (t in 1:nT) {
        for (s in 1:nEntangleStages){
          woundArrayList[[period]][i,,,s,t] <- woundProbMatrix(
            iTheta[i,], eps.i[,i,t],
            w.calf =  c(0,0,0,1,0)[s],
            resting = c(0,0,0,0,1)[s],
            age =     c(1,3,5,5,5)[s],
            reduce.iE = iE.change[[period]][s],
            reduce.iV = iV.change[[period]][s]*iV.change.t[[period]][s]^t)
        }}}}
  
  woundArray <- abind(woundArrayList,along=1)
}
