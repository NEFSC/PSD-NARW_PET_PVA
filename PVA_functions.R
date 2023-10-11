#object classes and functions for NARW PVA


#------------------------------------------------------------------------------#
# Individual whale class----
#------------------------------------------------------------------------------#
# a constructor function for the "whale" class
whale <- function(stage, wound = "fine", z = "alive") {
  # integrity checks 
  if (!(stage %in% stages)) {
    if (stage %in% 1:nStages)
      stage <- stages[stage]
    else
      stop("stage must be valid whale stage label (", paste(stages, collapse = ", "), ") or integer between 1 and ", nStages)
  }
  if (!(wound %in% woundStates))
    stop("wound must be valid severe wound state (", paste(woundStates, collapse = ", "), ")")
  if (!(z %in% liveDeadStates))
    stop("z must be in", liveDeadStates)
  value <- list(stage = factor(stage, levels = stages), 
                alive = factor(z, levels = liveDeadStates), 
                wound = factor(wound, levels = woundStates))
  attr(value, "class") <- "whale"
  value
}

# Determines whether this whale survives this year, and if it dies, of what
survive.whale <- function(obj, alpha, eps.m, i.ref.yr = 2019) {

  S <- switch (as.character(wound.whale(obj)),
               
               fine = exp(-exp(log(alpha["Mu.mO"]))),
               
               ent = exp(-exp(
                 log(alpha[["Mu.mE"]]) +
                   alpha[["a.mE.age"]]*(age.whale(obj = obj)-5) +
                   alpha[["a.mE.calf"]]*repro.whale(obj = obj) +
                   alpha[["a.mE.rest"]]*rest.whale(obj = obj) 
                   # alpha[["a.mE.regime"]]*as.numeric(i.ref.yr > 2010) + 
                   #alpha[["a.mE.regime2"]]*as.numeric(i.ref.yr > 2013)
                   # eps.m[["E"]]
                 )),
               
               vessel = exp(-exp(
                 log(alpha[["Mu.mV"]]) +
                   alpha[["a.mV.age"]]*(age.whale(obj = obj)-5) +
                   alpha[["a.mV.calf"]]*repro.whale(obj = obj) +
                   alpha[["a.mV.rest"]]*rest.whale(obj = obj) 
                   #alpha[["a.mV.regime"]]*as.numeric(i.ref.yr > 2010) +
                   #alpha[["a.mV.regime2"]]*as.numeric(i.ref.yr > 2013)
                   # eps.m[["V"]]
                 ))
               )
  lives <- rbinom(1, 1, S)
  if (lives){
    state <- 1
  } else {
    state <- switch (as.character(wound.whale(obj)),
                     fine = 4,
                     ent = 2,
                     vessel = 3)
  }
  liveDeadStates[state]
}

# Updates severe wound state for this whale and this year
updateWound.whale <- function(obj, woundProb) {
  sample(woundStates, 1, 
         prob = woundProb[, wound.whale(obj, how = 'short'), 
                          stagesToES[stage.whale(obj, "integer")]])
}

# Updates stage for this whale and this year
updateStage.whale <- function(obj, beta, eps.r, prey, N, kappa, B.ref.yr = 2019) {
  B <- plogis(beta[1:nReproStages] + 
                beta["b.prey1"] * prey[1] + beta["b.prey2"] * prey[2] + 
                as.numeric(beta["b.regime2"]) * as.numeric(B.ref.yr >= 2013) %*%
                   matrix(c(1,1,1,1,1,1,0),ncol = nReproStages,nrow = 1) +
                #beta["b.ent"] * (entangled.whale(obj)) + 
                #beta["b.ves"] * (struck.whale(obj)) +
                beta["b.inj"] * (entangled.whale(obj)) + 
                beta["b.inj"] * (struck.whale(obj)) +
                eps.r)
  names(B) <- c(5:10, "W")
  mat <- stageProbMatrix(B, kappa)
  sample(stages, 1, prob = mat[ , stage.whale(obj)])
}

# Transition whale's state from time t to t+1
update.whale <- function(obj, alpha, eps.m, beta, eps.r, prey, N, kappa, 
                         woundProb, B.ref.yr = 2019) {
  # If whale is already dead, don't change its state
  if (!alive.whale(obj))
    return(obj)
  out <- obj
  out$wound <- factor(updateWound.whale(obj, woundProb), levels = woundStates)
  # need the survival to be of the new wound updated "out" 
  out$alive <- factor(survive.whale(out, alpha, eps.m), 
                      levels = liveDeadStates)
  # If whale is now dead, don't change its state further
  if (!alive.whale(out))
    return(out)
  out$stage <- factor(updateStage.whale(obj, beta, eps.r, prey, N, kappa, B.ref.yr),
                      levels = stages) 
  return(out)
}

# Checks if this whale is a juvenile (TRUE/FALSE) 
juv.whale <- function(obj) {
  (obj$stage %in% c('F1', 'F2', 'F3', 'F4', 
                    'M1', 'M2', 'M3', 'M4'))
}

# Checks if this whale is female (TRUE/FALSE) 
fem.whale <- function(obj) {
  startsWith(as.character(obj$stage), 'F')
}

# Function for determining if this whale is alive
alive.whale <- function(obj) {
  obj$alive=="alive"
}

# Function for determining why this whale died
cause.whale <- function(obj, how = c("short", "long")) {
  how <- match.arg(how)
  switch(how,
         short = obj$alive,
         long = liveDeadStatesLong[which(liveDeadStates==obj$alive)])
}

# Function for determining if this whale just gave birth
repro.whale <- function(obj) {
  ifelse (obj$stage=='FC', 1, 0) # or however we're representing this stage
}

# Function for determining if this whale just weaned a calf
rest.whale <- function(obj) {
  ifelse (obj$stage=='FR', 1, 0) # or however we're representing this stage
}

# Function for determining this whale's stage (either printed short, long, or just position)
stage.whale <- function(obj, how = c('short', 'long', 'integer')) {
  how <- match.arg(how)
  switch(how,
         short = as.character(obj$stage),
         long = stagesLong[which(stages==obj$stage)],
         integer = which(stages==obj$stage))
}

# Function for determining this whale's survival age 
age.whale <- function(obj) {
  age1 <- stage.whale(obj = obj, "integer")
  if (fem.whale(obj = obj)) {
    if (age1 > 5)
      return(5)
    else
      return(age1)
  } 
  # males start at "M1" (stages==14) but need to subtract
  # number of female ages/stages to get the math right
  else 
    return(age1 - nFemStages)
}

# Function for determining this whale's reproductive stage for the purpose of 
# calculating survival rate
surv.stage.whale <- function(obj) {
  stage1 <- stage.whale(obj = obj, how = "short")
  age1 <- stage.whale(obj = obj, how = "integer")
  if (fem.whale(obj = obj)) {
    if (age1 < 5)
      return("alphaBreed1")
    return(switch(stage1,
                  FC = "alphaBreed3",
                  FR = "alphaBreed4",
                  FW = "alphaBreed5",
                  "alphaBreed2"))
  } 
  else 
    return("alphaBreed1")
}

# Function for determining this whale's wound state (either printed short or long)
wound.whale <- function(obj, how = c('short', 'long')) {
  how <- match.arg(how)
  switch(how,
         short = obj$wound,
         long = woundStatesLong[which(woundStates==obj$wound)])
}

# Function for determining this whale is severely entangled
entangled.whale <- function(obj) {
  ifelse(obj$wound=="ent", 1, 0)
}

# Function for determining this whale is severely wounded from vessel strike
struck.whale <- function(obj) {
  ifelse(obj$wound=="vessel", 1, 0)
}

# Function for determining this whale is severely wounded from other
other.whale <- function(obj) {
  ifelse(obj$wound=="other", 1, 0)
}

# Output of this whale's current state
print.whale <- function(obj, how = c('short', 'long')) {
  how <- match.arg(how)
  cat(cause.whale(obj, how = how), "\n")
  cat("Stage:", stage.whale(obj, how = how), "\n")
  cat("Wound State:", wound.whale(obj, how = how), "\n")
}



#------------------------------------------------------------------------------#
# Functions for list of whales----
#------------------------------------------------------------------------------#
nAlive <- function(whales) {
  sum(sapply(whales, alive))
}

nDead <- function(whales) {
  sum(sapply(whales, function(x) 1 - alive(x)))
}

nDeadCause <- function(whales) {
  table(sapply(whales, function(x) x$alive), exclude = "alive")
}

nRepro <- function(whales) {
  sum(sapply(whales, repro))
}

nEntangled <- function(whales) {
  sum(sapply(whales, entangled))
}

nStruck <- function(whales) {
  sum(sapply(whales, struck))
}

nOther <- function(whales) {
  sum(sapply(whales, other))
}

nFem <- function(whales) {
  sum(sapply(whales, fem))
}

nJuv <- function(whales) {
  sum(sapply(whales, juv.whale))
}

nAdult <- function(whales) {
  sum(sapply(whales, function(x) 1 - calf(x)))
}

stageStructure <- function(whales) {
  as.matrix(table(factor(sapply(whales, stage, how = 'short'), levels = stages)))
}

stageDeathStructure <- function(whales) {
  as.matrix(table(factor(sapply(whales, stage, how = 'short'), levels = stages),
                  factor(sapply(whales, cause, how = 'short'), levels = deadStates)))
}

printWhales <- function(whales, how = c('short', 'long')) {
  how <- match.arg(how)
  for (i in 1:length(whales)){
    print(whales[[i]], how)
    cat('\n')
  }
}
#------------------------------------------------------------------------------#
# Other/general functions----
#------------------------------------------------------------------------------#
# Creates matrix of wound state transition probabilities
woundProbMatrix <- function(iTheta, # injury parameters from COD model
                            eps.i, # temporal random effects (iE, iV)
                            w.calf = 0, # female with calf
                            resting = 0, # female resting from calf raising
                            age = 5, #adult (can range from 1-5)
                            # reduce the rate to some proportion of baseline
                            reduce.iV = 1, reduce.iE = 1,
                            # reference time period from the retrospective estimates
                            # only relevant for pre/post regime years (2010 or 2013)
                            i.ref.yr = 2019) {
  regime1 <- as.numeric(i.ref.yr > 2010)
  regime2 <- as.numeric(i.ref.yr > 2013)
  
  mat <- matrix(0, nWoundStates, nWoundStates, 
                dimnames = list(woundStates, woundStates))
  
  iV <- exp(log(iTheta[,"Mu.iV"]) + iTheta[,"a.iV.age"]*(age-5) + iTheta[,"a.iV.calf"]*w.calf +
              iTheta[,"a.iV.rest"]*resting +
              iTheta[,"a.iV.regime2"]*regime2 + eps.i[["V"]]) * reduce.iV
  
  iE <- exp(log(iTheta[,"Mu.iE"]) + iTheta[,"a.iE.age"]*(age-5) + iTheta[,"a.iE.calf"]*w.calf +
              iTheta[,"a.iE.rest"]*resting +
              iTheta[,"a.iE.regime2"]*regime2 + eps.i[["E"]]) * reduce.iE
  
  pinj <- 1-exp(-(iV + iE))
  pinj_E <- iE / (iV + iE)
  
  # col = t, row = t+1
  mat[1, 1] <- 1 - pinj
  mat[2, 1] <- pinj*pinj_E
  mat[3, 1] <- pinj*(1-pinj_E)
  mat[1, 2] <- 1 - iTheta[,"psiE"]
  mat[2, 2] <- iTheta[,"psiE"]
  mat[1, 3] <- 1
  mat
}

# Creates matrix of stage transition probabilities (given survival)
stageProbMatrix <- function(Bs, kappa) {
  names(Bs) <- reproStages
  mat <- matrix(0, nStages, nStages, dimnames = list(stages, stages))
  mat["F2", "F1"] <- mat["F3", "F2"] <- mat["F4", "F3"] <- 1
  mat["F5", "F4"] <- 1
  mat["F6", "F5"] <- 1 - Bs["5"]
  mat["F7", "F6"] <- 1 - Bs["6"]
  mat["F8", "F7"] <- 1 - Bs["7"]
  mat["F9", "F8"] <- 1 - Bs["8"]
  mat["F10+", "F9"] <- 1 - Bs["9"]
  mat["F10+", "F10+"] <- 1 - Bs["10"]
  mat["FC", "F5"] <- Bs["5"] * (1 - kappa)
  mat["FC", "F6"] <- Bs["6"] * (1 - kappa)
  mat["FC", "F7"] <- Bs["7"] * (1 - kappa)
  mat["FC", "F8"] <- Bs["8"] * (1 - kappa)
  mat["FC", "F9"] <- Bs["9"] * (1 - kappa)
  mat["FC", "F10+"] <- Bs["10"] * (1 - kappa)
  mat["FR", "F5"] <- Bs["5"] * kappa
  mat["FR", "F6"] <- Bs["6"] * kappa
  mat["FR", "F7"] <- Bs["7"] * kappa
  mat["FR", "F8"] <- Bs["8"] * kappa
  mat["FR", "F9"] <- Bs["9"] * kappa
  mat["FR", "F10+"] <- Bs["10"] * kappa
  mat["FR", "FC"] <- mat["FW", "FR"] <- 1
  mat["FW", "FW"] <- 1 - Bs["W"]
  mat["FC", "FW"] <- Bs["W"] * (1 - kappa)
  mat["FR", "FW"] <- Bs["W"] * kappa
  mat["M2", "M1"] <- mat["M3", "M2"] <- mat["M4", "M3"] <- 1
  mat["M5+", "M4"] <- mat["M5+", "M5+"] <- 1
  mat
}

# Creates matrix of stage transition probabilities (includes removing through 
# death and adding through birth)
stageProbMatrix2 <- function(Bs, kappa, S) {
  names(Bs) <- reproStages
  mat <- matrix(0, nStages, nStages, dimnames = list(stages, stages))
  mat["F2", "F1"] <- S[1] 
  mat["F3", "F2"] <- S[2] 
  mat["F4", "F3"] <- S[3]
  mat["F5", "F4"] <- S[4]
  mat["F6", "F5"] <- (1 - Bs["5"]) * S[5]
  mat["F7", "F6"] <- (1 - Bs["6"]) * S[5]
  mat["F8", "F7"] <- (1 - Bs["7"]) * S[5]
  mat["F9", "F8"] <- (1 - Bs["8"]) * S[5]
  mat["F10+", "F9"] <- (1 - Bs["9"]) * S[5]
  mat["F10+", "F10+"] <- (1 - Bs["10"]) * S[5]
  mat["FC", "F5"] <- Bs["5"] * (1 - kappa) * S[5]
  mat["FC", "F6"] <- Bs["6"] * (1 - kappa) * S[5]
  mat["FC", "F7"] <- Bs["7"] * (1 - kappa) * S[5]
  mat["FC", "F8"] <- Bs["8"] * (1 - kappa) * S[5]
  mat["FC", "F9"] <- Bs["9"] * (1 - kappa) * S[5]
  mat["FC", "F10+"] <- Bs["10"] * (1 - kappa) * S[5]
  mat["F1", "F5"] <- mat["M1", "F5"] <- Bs["5"] * (1 - kappa) * S[5] * 0.5
  mat["F1", "F6"] <- mat["M1", "F6"] <- Bs["6"] * (1 - kappa) * S[5] * 0.5
  mat["F1", "F7"] <- mat["M1", "F7"] <- Bs["7"] * (1 - kappa) * S[5] * 0.5
  mat["F1", "F8"] <- mat["M1", "F8"] <- Bs["8"] * (1 - kappa) * S[5] * 0.5
  mat["F1", "F9"] <- mat["M1", "F9"] <- Bs["9"] * (1 - kappa) * S[5] * 0.5
  mat["F1", "F10+"] <- mat["M1", "F10+"] <- Bs["10"] * (1 - kappa) * S[5] * 0.5
  mat["FR", "F5"] <- Bs["5"] * kappa * S[5]
  mat["FR", "F6"] <- Bs["6"] * kappa * S[5]
  mat["FR", "F7"] <- Bs["7"] * kappa * S[5]
  mat["FR", "F8"] <- Bs["8"] * kappa * S[5]
  mat["FR", "F9"] <- Bs["9"] * kappa * S[5]
  mat["FR", "F10+"] <- Bs["10"] * kappa * S[5]
  mat["FR", "FC"] <- S[6] 
  mat["FW", "FR"] <- S[7]
  mat["FW", "FW"] <- (1 - Bs["W"]) * S[8]
  mat["FC", "FW"] <- Bs["W"] * (1 - kappa) * S[8]
  mat["F1", "FW"] <- mat["M1", "FW"] <- Bs["W"] * (1 - kappa) * S[8] * 0.5
  mat["FR", "FW"] <- Bs["W"] * kappa * S[8]
  mat["M2", "M1"] <- S[1] 
  mat["M3", "M2"] <- S[2] 
  mat["M4", "M3"] <- S[3]
  mat["M5+", "M4"] <- S[4] 
  mat["M5+", "M5+"] <- S[9]
  mat
}

# Bias-corrected confidence intervals
bcCI <- function(x, probs, pointEst = mean(x, na.rm = T)) {
  z0 <- qnorm(sum((x)<pointEst, na.rm = T)/length(x))
  bcCI <- quantile(x, pnorm(2*z0+qnorm(probs)), na.rm=TRUE)
  return(bcCI)
}

# Expected minimum population size, given Ntot
emp <- function(Ntot) {
  emp <- array(NA, c(nBoot, nT, nS))
  for (t in 1:nT) {
    emp[, t, ] <- apply(Ntot[ , , 1:t, , drop=FALSE], c(1, 4), min)
  }
  return(emp)
}

# Compute the total number of whales in a subset of stages, such as adult
# females
NtotSubset <- function(N, stageSet = adFemStages) {
  apply(N[,,,,stageSet, drop = FALSE], 1:4, sum)
}


# For post-hoc computation of probability of quasi-extinction. Runs all
# scenarios. Ntot can be actual total or a precomputed subset (see previous
# function). Is there a way to run this without the quadruple-nested for loop? 
# Doesn't really matter, runs pretty fast.
compPQE <- function(Ntot, thresh = thresholds) {
  nThresholds <- length(thresh)
  PQE <- array(0, c(nBoot, nT, nS, nThresholds), 
               list(NULL, NULL, scenarios, thresh))
  for (sc in 1:nS)
    for (i in 1:nBoot) {
      for (j in 1:nRep) {
        for (th in 1:nThresholds) {
          t <- which(Ntot[i, j, , sc, drop = FALSE]<thresh[th])[1]
          if (!is.na(t))
            PQE[i, t:nT, sc, th] <- PQE[i, t:nT, sc, th] + 1/nRep
        }
      }
    }
  return(PQE)
}

compPDecline <- function(Ntot, fraction = c(0.7,0.5,0.2)){
  PDecline <- array(0, c(nBoot, nT, nS, length(fraction)), 
               list(NULL, NULL, scenarios, fraction))
  for (sc in 1:nS)
    for (i in 1:nBoot) {
      for (j in 1:nRep) {
        for (th in 1:length(fraction)) {
          t <- which(Ntot[i,j,,sc,drop=FALSE] < fraction[th]*Ntot[i,j,1,sc])[1]
          if (!is.na(t))
            PDecline[i, t:nT, sc, th] <- PDecline[i, t:nT, sc, th] + 1/nRep
        }
      }
    }
  return(PDecline)
}

compPIncrease <- function(Ntot, year=35, factor=2){
  PIncrease <- array(0, c(nBoot, nS), 
                    list(NULL, scenarios))
  for (sc in 1:nS)
    for (i in 1:nBoot) {
      for (j in 1:nRep) {
        PIncrease[i, sc] <- as.numeric(Ntot[i,j,year,sc] > factor*Ntot[i,j,1,sc])
      }
    }
  return(PIncrease)
}

nBootKeep <- function(posterior){
  if(exists("nBoot")){
    return(floor(seq(1,nrow(posterior),length.out=nBoot)))
  } else {
    print("nBoot missing, using 1,000 samples")
    nBoot <- 1000
    return(floor(seq(1,nrow(posterior),length.out=nBoot)))
  }
}

getCI <- function(x,prob=0.95,rd=4,na_rm=T,stat=c("med","mean")[1]){
  probs <- c(0.50,0+((1-prob)/2),1-((1-prob)/2))
  if (stat=="mean"){
    return(round(c(mean(x),quantile(x,probs=probs[-1],na.rm=na_rm)),rd))
  } else {
    return(round(quantile(x,probs=probs,na.rm=na_rm),rd))
  }
}

# Run PVA function, for a single scenario 
# Inputs: 
#   params - list with parameter inputs, currently N0 (nStages x nBoot matrix of 
#     starting pop sizes), betas (nBoot x 12 array of reproduction parameters), 
#     kappa (nBoot probability of early calf loss), alphas (nBoot x 9 array of 
#     survival parameters), prey (nBoot x nT array of prey values, 
#     nu (nBoot x nT matrix of survival random values), epsilon (nBoot x nT 
#     matrix of repro random values),  wound0 (list of wound values for 
#     each starting whale), woundProb (array of severe wound probabilities),
#     fractions (array of proportion of mortality due each cause), and
#     prop (array of proportional changes in CSM)
#   nBoot - number of parametric bootstrap runs to loop through
#   nRep - number of Monte Carlo replications to loop through
#   ceiling_N - maximum number of whales the population can hold (ceiling DD)
#   verbose - integer 1 - 5
#     1: no display during function run
#     2: display bootstrap runs
#     3: display bootstrap runs and Monte Carlo replications
#     4: display bootstrap runs, Monte Carlo replications, years, and population summary
#     5: display bootstrap runs, Monte Carlo replications, years, population summary, and individual whales
runPVA <- function(params, nBoot = 1000, nRep = 1, nT = 100, ceiling_N = 5000, 
                   verbose = 1) {
  # Data structures 
  N <- array(NA, c(nBoot, nRep, nT, nStages), dimnames = list(NULL, NULL, NULL, stages))
  Ntot <- array(NA, c(nBoot, nRep, nT))
  PQE <- array(0, c(nBoot, nT, nThresholds), dimnames = list(NULL, NULL, thresholds))
  Ndead <- array(0, c(nBoot, nRep, nT, nStages, nDeadStates),
                 list(NULL, NULL, NULL, stages, deadStates))
  Nborn <- array(0, c(nBoot, nRep, nT))
  propEntangled <- array(NA, c(nBoot, nRep, nT))
  propStruck <- array(NA, c(nBoot, nRep, nT))
  propOther <- array(NA, c(nBoot, nRep, nT))
  for (i in 1:nBoot) {
    if (verbose > 1)
      cat("Bootstrap Run:", i, "of", nBoot, 'at', format(Sys.time()), "\n")
    N0rep <- rep(stages, params$N0[i, ])
    wound0rep <- rep(woundStates, apply(params$wound0[i, , ], 2, sum))
    for (j in 1:nRep) {
      if (verbose > 2)
        cat("\tMonte Carlo Run:", j, "of", nRep, 'at', format(Sys.time()), "\n")
      # Need somehow to come up with starting conditions in year 1
      whales <- mapply(whale, stage = N0rep, wound = wound0rep, SIMPLIFY = FALSE)
      for (t in 1:nT) {
        if (verbose > 3) {
          cat("Year:", t, "of", nT, "\n")
          cat(length(whales), "whales\n")
        }
        if (length(whales) == 0) {
          N[i, j, t:nT, ] <- 0
          Ntot[i, j, t:nT] <- 0
          break
        }
        whales <- lapply(whales, update, alpha = params$alphas[i, ], 
                         #iTheta = params$iTheta[i, ],
                         beta = params$betas[i, ], 
                         woundProb = params$woundProb[i,,,,t],
                         #nu = params$nu[i, t], 
                         eps.m = params$eps.m[ , i, t],
                         eps.r = params$eps.r[i, t],
                         prey = params$prey[i, t, ],
                         N = length(whales),
                         kappa = params$kappa[i]
                         #fractions = params$fractions[i, , ],
                         #prop = params$prop[ , , t]
        )
        if (verbose > 4)
          printWhales(whales, "long")
        currentDead <- which(sapply(whales, function(x) !alive(x)))
        if (length(currentDead)>0){
          #print("Dead:")
          #printWhales(whales[currentDead])
          Ndead[i, j, t, , ] <- (stageDeathStructure(whales[currentDead]))
          #print(Ndead[i, j, t, , ])
          whales <- whales[-currentDead]
        }
        nInd <- length(whales)
        if (nInd == 0){
          N[i, j, t:nT, ] <- 0
          Ntot[i, j, t:nT] <- 0
          break
        }
        nRecruit <- nRepro(whales)
        Nborn[i, j, t] <- nRecruit
        if (nRecruit > 0) {
          fem <- rbinom(nRecruit, 1, 0.5)
          for (ind in 1:nRecruit)
            # may indicate stages differently
            whales[[ind + nInd]] <- whale(ifelse(fem[ind], 'F1', 'M1'), 
                                          sample(woundStates, 1, 
                                                 prob = params$woundProb[i,,1,1,t]))
          nInd <- nInd + nRecruit
        }
        # Ceiling DD, put in just in case
        if (nInd > ceiling_N) {
          whales <- whales[1:ceiling_N]
          nInd <- ceiling_N
          # to speed up projections with fast growing N
          # Ntot[i, j, t:nT] <- ceiling_N
          # break
        }
        N[i, j, t, ] <- as.vector(stageStructure(whales))
        Ntot[i, j, t] <- nInd
        propEntangled[i, j, t] <- nEntangled(whales) / nInd
        propStruck[i, j, t] <- nStruck(whales) / nInd
        propOther[i, j, t] <- nOther(whales) / nInd
      }
      for (th in 1:nThresholds) {
        t <- which(Ntot[i, j, ]<thresholds[th])[1]
        if (!is.na(t))
          PQE[i, t:nT, th] <- PQE[i, t:nT, th] + 1/nRep
      }
    } #j, rep
  } #i, boot
  return(list(PQE = PQE, N = N, Ntot = Ntot, Ndead = Ndead, Nborn = Nborn,
              propEntangled = propEntangled, propStruck = propStruck,
              propOther = propOther))
}

#------------------------------------------------------------------------------#
# Generic functions----
#------------------------------------------------------------------------------#
survive <- function(obj, ...) {
  UseMethod("survive")
}

alive <- function(obj, ...) {
  UseMethod("alive")
}

repro <- function(obj, ...) {
  UseMethod("repro")
}

stage <- function(obj, ...) {
  UseMethod("stage")
}

ent <- function(obj, ...) {
  UseMethod("ent")
}

entangled <- function(obj, ...) {
  UseMethod("entangled")
}

struck <- function(obj, ...) {
  UseMethod("struck")
}

other <- function(obj, ...) {
  UseMethod("other")
}

fem <- function(obj, ...) {
  UseMethod("fem")
}

juv <- function(obj, ...) {
  UseMethod("juv")
}

cause <- function(obj, ...) {
  UseMethod("cause")
}

#------------------------------------------------------------------------------#
# Parallel run----
#------------------------------------------------------------------------------#

runPVA.par <- function(params, nRep = 1, nT = 100, ceiling_N = 5000, 
                   verbose = 1) {
  # Data structures 
  N <- array(NA, c(nRep, nT, nStages), dimnames = list(NULL, NULL, stages))
  Ntot <- array(NA, c( nRep, nT))
  PQE <- array(0, c( nT, nThresholds), dimnames = list(NULL, thresholds))
  Ndead <- array(0, c(nRep, nT, nStages, nDeadStates),
                 list(NULL, NULL, stages, deadStates))
  Nborn <- array(0, c(nRep, nT))
  propEntangled <- array(NA, c(nRep, nT))
  propStruck <- array(NA, c(nRep, nT))
  propOther <- array(NA, c(nRep, nT))
  
  # could use source("PVA_functions.R")
  # and clusterExport(cl, c("params","nBoot","nRep","nT","ceiling_N","verbose"))
  
  # for (i in 1:nBoot) {
  #   if (verbose > 1)
  #     cat("Bootstrap Run:", i, "of", nBoot, 'at', format(Sys.time()), "\n")
    N0rep <- rep(stages, params$N0)
    wound0rep <- rep(woundStates, apply(params$wound0, 2, sum))
    for (j in 1:nRep) {
      if (verbose > 2)
        cat("\tMonte Carlo Run:", j, "of", nRep, 'at', format(Sys.time()), "\n")
      # Need somehow to come up with starting conditions in year 1
      whales <- mapply(whale, stage = N0rep, wound = wound0rep, SIMPLIFY = FALSE)
      for (t in 1:nT) {
        if (verbose > 3) {
          cat("Year:", t, "of", nT, "\n")
          cat(length(whales), "whales\n")
        }
        if (length(whales) == 0) {
          N[j, t:nT, ] <- 0
          Ntot[j, t:nT] <- 0
          break
        }
        whales <- lapply(whales, update, alpha = params$alphas, 
                         beta = params$betas, 
                         woundProb = params$woundProb[,,,t],
                         eps.m = params$eps.m[ , t],
                         eps.r = params$eps.r[t],
                         prey = params$prey[t, ],
                         N = length(whales),
                         kappa = params$kappa,
                         B.ref.yr = params$B.ref.yr
        )
        if (verbose > 4)
          printWhales(whales, "long")
        currentDead <- which(sapply(whales, function(x) !alive(x)))
        if (length(currentDead)>0){
          #print("Dead:")
          #printWhales(whales[currentDead])
          Ndead[j, t, , ] <- (stageDeathStructure(whales[currentDead]))
          #print(Ndead[i, j, t, , ])
          whales <- whales[-currentDead]
        }
        nInd <- length(whales)
        if (nInd == 0){
          N[j, t:nT, ] <- 0
          Ntot[j, t:nT] <- 0
          break
        }
        nRecruit <- nRepro(whales)
        Nborn[j, t] <- nRecruit
        if (nRecruit > 0) {
          fem <- rbinom(nRecruit, 1, 0.5)
          for (ind in 1:nRecruit)
            # may indicate stages differently
            whales[[ind + nInd]] <- whale(ifelse(fem[ind], 'F1', 'M1'), 
                                          sample(woundStates, 1, 
                                                 prob = params$woundProb[,1,1,t]))
          nInd <- nInd + nRecruit
        }
        # Ceiling DD, put in just in case
        if (nInd > ceiling_N) {
          whales <- whales[1:ceiling_N]
          nInd <- ceiling_N
        }
        N[j, t, ] <- as.vector(stageStructure(whales))
        Ntot[j, t] <- nInd
        propEntangled[j, t] <- nEntangled(whales) / nInd
        propStruck[j, t] <- nStruck(whales) / nInd
        propOther[j, t] <- nOther(whales) / nInd
      }
      
      PQE <- array(0, c(nT, nThresholds), dimnames = list(NULL, thresholds))
      for (th in 1:nThresholds) {
        t <- which(Ntot[j, ]<thresholds[th])[1]
        if (!is.na(t))
          PQE[t:nT, th] <- PQE[t:nT, th] + 1/nRep
      }
    } #j, rep
  #} #i, boot
  return(list(PQE = PQE, N = N, Ntot = Ntot, Ndead = Ndead, Nborn = Nborn,
              propEntangled = propEntangled, propStruck = propStruck,
              propOther = propOther))
}
