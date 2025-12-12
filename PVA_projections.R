#------------------------------------------------------------------------------#
# Population Evaluation Tool (PET) for the North Atlantic Right Whale (NARW)
# Version 1.0

# Structure of input files (actual names may vary):
#   PVA_population_t0.Rdata: nBoot or more samples of starting stage and wound
#      distributions.
#   NARW_posteriors_REPRO.csv: nBoot (or more) samples of parameters from 
#     integrated survival/reproduction model. 
#   NARW_posteriors_MORT.csv: nBoot (or more) samples of fractions of mortality and 
#     wounding probabilities.
#   NARW_food_covariates_1986-2019.Rdata: Historical data on prey availability, not yet 
#     sampled, averaged, or weighted.
#------------------------------------------------------------------------------#
library(zoo);library(abind)
out_drive <- "C:\\temp\\"

#------------------------------------------------------------------------------#
# Basic settings----
#------------------------------------------------------------------------------#
version <- "1.2"

nBoot <- 1000#, number of bootstrap runs
nRep  <- 1   #, number of replications (Monte Carlo loop). 
nT <- 100 # number of years

## Reproductive stages----
stages <- c('F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 
            'F9', 'F10+', 'FC', 'FR', 'FW', 
            'M1', 'M2', 'M3', 'M4', 'M5+')
stagesLong <- c('Female Prebreeder 0.5 Years Old', 'Female Prebreeder 1.5 Years Old', 
                'Female Prebreeder 2.5 Years Old', 'Female Prebreeder 3.5 Years Old', 
                'Female Prebreeder 4.5 Years Old', 'Female Prebreeder 5.5 Years Old', 
                'Female Prebreeder 6.5 Years Old', 'Female Prebreeder 7.5 Years Old', 
                'Female Prebreeder 8.5 Years Old', 'Female Prebreeder 9.5+ Years Old', 
                'Female With Calf', 'Female Resting', 'Female Breeder (No Calf)', 
                'Male 0.5 Years Old', 'Male 1.5 Years Old', 'Male 2.5 Years Old', 
                'Male 3.5 Years Old', 'Male 4.5+ Years Old')
nStages <- length(stages)
femStages <- stages[grep("F", stages)]
nFemStages <- length(femStages)
adFemStages <- femStages[-(1:4)]

### Reproductive stage effects structure on B and S----
reproStages <- c(5:10, "W")
nReproStages <- length(reproStages)
survAges <- 1:5
survStages <- 1:5
nSurvAges <- length(survAges)
nSurvStages <- length(survStages)

## Wound/death states----
woundStates <- c("fine", "ent", "vessel") #, "other")
woundStatesLong <- c("Not Wounded (Or Not Severely)", "Severe Entanglement Wound",
                     "Severe Vessel Wound") #, "Severe Other Wound")
nWoundStates <- length(woundStates)

deadStates <- c("ent", "vessel", "other")
deadStatesLong <- c("Recently Dead Entanglement", "Recently Dead Vessel",
                    "Recently Dead Other")
nDeadStates <- length(deadStates)

liveDeadStates <- c("alive", deadStates)
liveDeadStatesLong <- c("Alive", deadStatesLong)
nLiveDeadStates <- nDeadStates + 1

entangleStages <- c("Calf", "Juvenile", "Adult", "FC", "FR") 
nEntangleStages <- length(entangleStages)
stagesToES <- c(rep(entangleStages[-c(4:5)], c(1, 3, 9)), #Females
                rep(entangleStages[-c(4:5)], c(1, 3, 1))) #Males
# females currently/recently w/ calf are entanglement specific states
stagesToES[stages %in% c("FC","FR")] <- c("FC","FR")

## Extinction thresholds----
thresholds <- c(1, 10, 50, 100)
nThresholds <- length(thresholds)

ceiling_N <- 10000 # Just in case
verbose <- 2 #integer 0 to 5
seed <- 2024 #opportunity to replicate random sampling
set.seed(seed)

#------------------------------------------------------------------------------#
# Read in the functions----
#------------------------------------------------------------------------------#
source("PVA_functions.R")
source("PVA_scenarios.R")

#------------------------------------------------------------------------------#
# Read in parameters, set up starting values, etc.----
#------------------------------------------------------------------------------#

## posterior distributions from reproduction and mortality models----
post_repro <- read.csv("inputs/NARW_posteriors_REPRO.csv")
post_mort  <- read.csv("inputs/NARW_posteriors_MORT.csv")

## match Greek letter convention for mortality model
mort.coeff <- grep(paste(c("beta.m","beta.i"),collapse="|"),names(post_mort),value=T)
names(post_mort)[match(mort.coeff,names(post_mort))] <- gsub("beta","a",mort.coeff)
  
## reduce size to number of bootstrap runs per projection----
post_repro <- post_repro[nBootKeep(post_repro),]
post_mort  <-  post_mort[nBootKeep(post_mort),]

## natural mortality (set to 0)----
# based on Dureuil & Froese (2021) and t_max = 268 (bowhead whales)
post_mort$Mu.mO <- 0 #rgamma(nrow(post_mort),1,99)

## starting stage/state distributions----
load("inputs/PVA_population_t0.Rdata")
N0 <- N0[nBootKeep(N0),]
dimnames(N0) <- list(NULL,stages)
wound0 <- wound0[nBootKeep(wound0),,]
dimnames(wound0) <- list(NULL, stages, woundStates)

## prey index data----
load(file="inputs/NARW_food_covariates_1986-2019.Rdata")

## reproductive parameters----
post_repro$beta.ageW <- qlogis(post_repro$p.beta.ageW)
post_repro[, grep("beta.age.first",colnames(post_repro),value=TRUE)] <-
  post_repro[, grep("beta.age.first",colnames(post_repro),value=TRUE)] + post_repro[,"beta.ageW"]
betas <- post_repro[, c(grep("beta.age.first",colnames(post_repro),value=TRUE),
                        "beta.regime", "beta.regime.p2", 
                        "beta.p1", "beta.p2", 
                        "beta.inj"
                        )]
colnames(betas) <- c(5:10, "W", 
                     "b.regime","b.regime.prey2",
                     "b.prey1", "b.prey2", "b.inj") 

## calf loss parameter----
kappa <- post_repro[, "kappa"]

## survival/mortality parameters----
alphas <- post_mort[,c(
  "Mu.mO","Mu.mE", "Mu.mV",
  "a.mE.age", "a.mE.calf", "a.mE.rest",
  "a.mV.age", "a.mV.calf", "a.mV.rest"
)]

## injury rates/coefficients----
iTheta <- post_mort[,c(
  "Mu.iE","Mu.iV","psiE",
  "a.iE.age", "a.iE.calf", "a.iE.regime2", "a.iE.rest", 
  "a.iV.age", "a.iV.calf", "a.iV.regime2", "a.iV.rest" 
)]

## random effects----
## annual deviations in injury/mortality and reproduction
eps.i <- eps.m <- array(0, c(2, nBoot, nT),
                        list(c("E","V"), NULL, NULL))
eps.r <- eps.i[1,,]
eps <- list(random = list(eps.i,eps.m,eps.r),
            # [NOT USED] for using estimated epsilon values
            period = list(eps.i,eps.m,eps.r))
names(eps[[1]]) <- names(eps[[2]]) <- c("inj","mort","repro")
rm(list = c("eps.i","eps.m","eps.r"))
eps.period <- which(1991:2019 > 2013)

for (i in 1:nBoot) {
  # sample of year-specific deviations
  r.samp <- sample(1:length(eps.period),nT,replace=T)
  
  # injury
  eps[["random"]][["inj"]][1, i, ] <- rnorm(nT, 0, post_mort[i, "sigma.iE.t"])
  eps[["random"]][["inj"]][2, i, ] <- rnorm(nT, 0, post_mort[i, "sigma.iV.t"])

  # mortality
  eps[["random"]][["mort"]][1, i, ] <- 0 #rnorm(nT, 0, post_mort[i, "sigma.mE.t"])
  eps[["random"]][["mort"]][2, i, ] <- 0 #rnorm(nT, 0, post_mort[i, "sigma.mV.t"])

  # reproduction
  eps[["random"]][["repro"]][i, ] <- rnorm(nT, 0, post_repro[i, "sigma.B.t"])
}
eps.type <- "random"

#------------------------------------------------------------------------------#
# Scenario building----
#------------------------------------------------------------------------------#
load.prey <- FALSE
preyChange_v <- c(0.7,0.8,0.9,1.1,1.2,1.3)

if (!load.prey){
preyArrays <- list()

# post-2010 prey availability
preyArrays[["prey post2010"]] <- build.scenario.prey(
  plankton_w_food,
  ref_yrs = list(2010:2019),
  proj_yrs = list(1:nT)
  )
# historic fluctuations
preyArrays[["prey historic"]] <- build.scenario.prey(
  plankton_w_food,
  # good vs. bad decade
  #ref_yrs = rep(list(2001:2009,1990:2000),nT/20),
  ref_yrs = list(1990:2009)
  # decadal cycling
  # proj_yrs = split(1:nT,ceiling((1:nT)/10))
)

# add prey availability due to noise
for (a in 1:6){
  preyArrays[[a+2]] <- build.scenario.prey(
    plankton_w_food,
    ref_yrs = list(2010:2019),
    proj_yrs = list(1:nT), preyChange = preyChange_v[a]
    
  )
}
names(preyArrays)[(1:length(preyChange_v))+2] <- paste0(
  "prey post2010 ", preyChange_v*100,"%")
save(preyArrays,preyChange_v,file="./inputs/preyArrays_ALLscenarios.Rdata")
}
load(file="./inputs/preyArrays_ALLscenarios.Rdata")


woundArrays <- list()

#ent.reduce <- c(0,25)
ent.reduce <- c(seq(0,100,by=10),25)
# total entanglement reductions
start <- Sys.time()
for (e in ent.reduce){
  woundArrays[[paste0("reduce entanglement ",e,
                      "%; vessel strike up 0%; speed restrict 0%")]] <- build.scenario.wound(
    iTheta,
    eps.i = eps[[eps.type]][["inj"]],
    iE.change = list(rep((100-e)/100,5))
  )
}
(end <- Sys.time()-start)

# weak rope coverage (entanglement reductions)
weak.e <- c(50)
for (e in weak.e){
  woundArrays[[paste0("weak rope coverage ",e,
                      "%; vessel strike up 0%; speed restrict 0%")]] <- build.scenario.wound(
    iTheta,
    eps.i = eps[[eps.type]][["inj"]],
    iE.change = list(c(rep((100-(e*.9))/100,2),rep((100-(e*.6))/100,3)))
  )
}

# all vessel scenarios (regarding speed restrictions/temporal change)
ent.reduce <- c(0,25)
vess_t_change <- c(0.7,0,-0.3)
vess_speed <- c(1,0.75,0)

#entanglement reduction
for (e in ent.reduce){        
  #speed restrict
  for (v2 in vess_speed[list(1:3,1:2)[[match(e,ent.reduce)]]]){
  #traffic change
  for (v1 in vess_t_change[list(1:3,1:3,2)[[match(v2,vess_speed)]]]){  
    woundArrays[[paste0(
      "reduce entanglement ",e,"%; vessel strike up ",v1,
      "%; speed restrict ",(1-v2)*100,"%")]] <- build.scenario.wound(
      iTheta,
      eps.i = eps[[eps.type]][["inj"]],
      iE.change = list(rep(1-(e/100),5)), #iE = 0.50 of status quo rate
      iV.change = list(rep(v2,5)),
      iV.change.t = list(rep(1+(v1/100),5))
    )
  }
  }
}

#-----------------------------------#
# count scenarios and name them
scenario.dat <- data.frame(
  wound = names(woundArrays),
  prey = names(preyArrays)[1],
  B.regime = "repro post2010",
  ceiling_N = ceiling_N
)
#-----------------------------------#

# add historical prey scenarios
scenario.dat <- rbind(
  scenario.dat,
  data.frame(
    wound = names(woundArrays)[c(1,19)], #baseline1 & 100% vess reduce
    prey = names(preyArrays)[c(2)],
    B.regime = c("repro post2010"),
    ceiling_N = ceiling_N
  ))

# add prey accessibility scenarios
scenario.dat <- rbind(
  scenario.dat,
  data.frame(
    wound = names(woundArrays)[c(1)], #0% ent reduction
    prey = names(preyArrays)[(1:length(preyChange_v))+2],
    B.regime = c("repro post2010"),
    ceiling_N = ceiling_N
  ))


scenario.dat$names <- apply(scenario.dat[,c("wound","prey")],1,function(x) paste(x,collapse="; "))
scenarios <- scenario.dat$names
nS <- nrow(scenario.dat)

#------------------------------------------------------------------------------#
# Output data structures----
#------------------------------------------------------------------------------#
N <- array(NA, c(nBoot, nRep, nT, nS, nStages), 
           dimnames = list(NULL, NULL, NULL, scenarios, stages))
Ntot <- array(NA, c(nBoot, nRep, nT, nS))
PQE <- array(0, c(nBoot, nT, nS, nThresholds), 
             dimnames = list(NULL, NULL, scenarios, thresholds))
Ndead <- array(0, c(nBoot, nRep, nT, nS, nStages, nDeadStates),
               list(NULL, NULL, 1:nT, scenarios, stages, deadStates))
Nborn <- array(0, c(nBoot, nRep, nT, nS),
               list(NULL, NULL, 1:nT, scenarios))
propEntangled <- array(NA, c(nBoot, nRep, nT, nS))
propStruck <- array(NA, c(nBoot, nRep, nT, nS))
propOther <- array(NA, c(nBoot, nRep, nT, nS))


#------------------------------------------------------------------------------#
# Run simulation----
#------------------------------------------------------------------------------#
library(parallel)
library(foreach)
library(doParallel)
n.cores <- 7

start.scenario <- 1

for (scenario in start.scenario:nS) {
  if (verbose > 0)
    cat("Running", scenarios[scenario], "scenario\n")
    print(Sys.time())
  
  #new cluster processing
  cl<-makeCluster(n.cores)
  clusterExport(cl,varlist=ls())
  registerDoParallel(cl)
  
  print(system.time(
    temp <- foreach(i=1:nBoot) %dopar%
      runPVA.par(
        params= list(
          N0 = N0[i,], betas = as.matrix(betas)[i,], 
          eps.i = eps[[eps.type]][["inj"]][,i,], 
          eps.m = eps[[eps.type]][["mort"]][,i,], 
          eps.r = eps[[eps.type]][["repro"]][i,],
          alphas = as.matrix(alphas)[i,],
          woundProb = woundArrays[[scenario.dat$wound[scenario]]][i,,,,],
          prey = preyArrays[[scenario.dat$prey[scenario]]][i,,],
          B.ref.yr = c(2019,1999)[(scenario.dat$B.regime[scenario]=="repro historic") + 1],
          wound0 = wound0[i,,], 
          kappa = kappa[i]
          ),
        ceiling_N = scenario.dat$ceiling_N[scenario],
        nT = nT,
        nRep = nRep
        )))
  
  stopCluster(cl)
  gc()
  
  
  PQE[ , , scenario, ] <- abind(lapply(1:nBoot,function(x) temp[[x]]$PQE),along=0)
  N[ , , , scenario, ] <- abind(lapply(1:nBoot,function(x) temp[[x]]$N),along=0)
  Ntot[ , , , scenario] <- abind(lapply(1:nBoot,function(x) temp[[x]]$Ntot),along=0)
  Ndead[ , , , scenario, , ] <- abind(lapply(1:nBoot,function(x) temp[[x]]$Ndead),along=0)
  Nborn[ , , , scenario] <- abind(lapply(1:nBoot,function(x) temp[[x]]$Nborn),along=0)
  propEntangled[ , , , scenario] <- abind(lapply(1:nBoot,function(x) temp[[x]]$propEntangled),along=0)
  propStruck[ , , , scenario] <- abind(lapply(1:nBoot,function(x) temp[[x]]$propStruck),along=0)
  propOther[ , , , scenario] <- abind(lapply(1:nBoot,function(x) temp[[x]]$propOther),along=0)
  #save.image(paste0(out_drive,'pvaResults', version, letters[1], '.RData'))
}
