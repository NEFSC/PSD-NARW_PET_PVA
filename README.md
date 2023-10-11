# A Management-focused Population Viability Analysis (PVA) for North Atlantic Right Whales

This repository contains the data and scripts developed by the Population Evaluation Team (PET) subgroup [^readme-1] to assess the effect of potential management measures on the future population trajectory of North Atlantic right whales. The PET PVA consists of a core stage-structured model and several submodels (mortality, reproduction, entanglement, vessel strike, prey availability, and prey accessibility). The primary scripts are as follows:

[^readme-1]: <https://www.fisheries.noaa.gov/new-england-mid-atlantic/endangered-species-conservation/north-atlantic-right-whale-recovery-plan-northeast-us-implementation-team>

-   `PVA_functions.R`
-   `PVA_scenarios.R`
-   `PVA_projections.R`

The first script uses object-oriented programming to define functions and attributes consistent with the life-history of NARWs. The second script is used to design scenarios related to wounding rates and prey dynamics. The third script sets up the parameters of the population projection and simulates an `nBoot` number of replicated population trajectories that incorporate stochasticity and uncertainty.

## Inputs

The inputs folder contains various objects representing both data and posterior distributions from other models that were used to provide estimates of demographic parameters. These include:

-   `NARW_food_covariates_1986-2019.Rdata` = indices of prey availability
-   `NARW_posteriors_MORT.csv` = parameters related to injury/mortality
-   `NARW_posteriors_REPRO.csv` = parameters related to reproduction
-   `PVA_population_t0.Rdata` = starting population size in 2019

## Projections

Conditional on the coded structure of the population and the specified collection of parameters defining demographic attributes, simulated populations are projected through time (e.g., 100 years) to examine distributions of future population sizes and probabilities of extinction (or quasi-extinction). These outputs can be used to examine anticipated outcomes of management measures that are aimed at reducing threats to the population.

*This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.*
