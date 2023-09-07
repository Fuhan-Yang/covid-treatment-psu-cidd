###############################################################################
### 1. Preliminaries
###############################################################################
rm(list = ls())
### set up working dir to be github repo 
setwd('~/github/covid-treatment-psu-cidd/')

### necessary packages
## library for spline-based expansions
library(fda)
## package for iSplines:
library(splines2)
## package for plotting
library(reshape2)

### load functions
source("R/functions/traj.from.params-v7-TT.R")
source("R/functions/main_traj.R")
source("R/functions/run_simulation_functions.R")
source("R/functions/generate_beta_func2.R")
source('R/functions/US_scale_up.R')

#### run simulation.R and baseline_params.func should be under the same dir ###

### load the tv-flu-vac-cov ###
load('data/us_flu_vac/relative_flu_vac_perc_by_age.Rdata')

############ load initial parameters: params fixed before Mar 1 2023 ###########
### load params under three scenarios ###
load('R/params//initial_params_365.Rdata')
# load('R/params/initial_params_548.Rdata')
# load('R/params/initial_params_730.Rdata')

all_params = all_params_365

### needed metrics ###
primary_metrics = c('cumu_cases','cumu_hosp','cumu_death')

##### run multiple params as parameter combinations ########
treat_cover = c(0.137,seq(0,1,0.05))
vac_cover = c(0.49,seq(0,1,0.05))

all_param_combs_list = list()

for(j in 1:length(all_params)) {
  
  param_combs = expand.grid(treat = treat_cover,
                            vac = vac_cover,
                            cumu_cases = NA,
                            cumu_hosp = NA,
                            cumu_death = NA)
  
  for(i in 1:nrow(param_combs)) {
    
    treat_by_age = rep(param_combs$treat[i],8)
    names(treat_by_age) = c("tv-treat-cov-10_2","tv-treat-cov-20_2",
                            "tv-treat-cov-30_2","tv-treat-cov-40_2",
                            "tv-treat-cov-50_2","tv-treat-cov-60_2",
                            "tv-treat-cov-70_2","tv-treat-cov-80_2")
    
    summary_traj = run_single_param_comb(change_params = treat_by_age,
                                         change_vac_cover = param_combs$vac[i], 
                                         pop_by_age = get_pop_by_age(loc = 'RI'),
                                         vac_tv_type = 'scenarios',
                                         vac_cov_ts = flu_vac_by_age,
                                         summary_type = 'diff',
                                         sim_startday = as.numeric(as.Date('2025-03-01') - as.Date('2020-01-01')) + 1,
                                         baseline_params = all_params[[j]],
                                         scale_us = T,
                                         scaled_metrics = primary_metrics)
    
    param_combs$cumu_cases[i] = summary_traj['cumu_cases']
    param_combs$cumu_hosp[i] = summary_traj['cumu_hosp']
    param_combs$cumu_death[i] = summary_traj['cumu_death']
  
  }
  all_param_combs_list[[j]] = param_combs
  
}

names(all_param_combs_list) = names(all_params)
save(all_param_combs_list,file = 'R/Rdata/vac_treat_summaries_sim365_20250301.Rdata')

