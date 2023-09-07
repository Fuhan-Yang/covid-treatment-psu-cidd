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

### load params under three scenarios ###
# load('R/params/initial_params_548.Rdata')

### load the tv-flu-vac-cov ###
load('data/us_flu_vac/relative_flu_vac_perc_by_age.Rdata')

baseline_params = all_params_548$neu

#### load the baseline treat cov and vac cov for each state ####
load('R/Rdata/treat_cov_baseline_by_state.Rdata')
load('R/Rdata/vac_cov_baseline_by_state.Rdata')
data(state)

### only include 50 states ###
pct_vacc_by_state = 
  pct_vacc_by_state %>%
  filter(Location %in% state.abb) %>%
  rename(state = Location)
treat_cov_by_state_df = treat_cov_by_state_df[treat_cov_by_state_df$state %in% state.abb,]

treat_vac_state_df = plyr::join(data.frame(state = pct_vacc_by_state$state,
                                     opt_vac = pct_vacc_by_state$opt_pct_vacc),
                          data.frame(state = treat_cov_by_state_df$state,
                                     treat = treat_cov_by_state_df$cov/100),
                          by = 'state')

#### manually add national baseline: treat cov = 0.137, vac cov = 0.49 ####
treat_vac_state_df = rbind(treat_vac_state_df,
                           data.frame(state = 'US',
                                      opt_vac = 0.49,
                                      treat = 0.137))


####### The situations in next ten years for each state if do nothing ########
treat_vac_state_df$cumu_cases = NA
treat_vac_state_df$cumu_hosp = NA
treat_vac_state_df$cumu_death = NA

primary_metrics = c('cumu_cases','cumu_hosp','cumu_death')

for(i in 1:nrow(treat_vac_state_df)) {
  
  treat_by_age = rep(treat_vac_state_df$treat[i],8)
  names(treat_by_age) = c("tv-treat-cov-10_2","tv-treat-cov-20_2",
                          "tv-treat-cov-30_2","tv-treat-cov-40_2",
                          "tv-treat-cov-50_2","tv-treat-cov-60_2",
                          "tv-treat-cov-70_2","tv-treat-cov-80_2")
  
  summary_traj = run_single_param_comb(change_params = treat_by_age,
                                       change_vac_cover = treat_vac_state_df$opt_vac[i], 
                                       pop_by_age = get_pop_by_age(loc = 'RI'),
                                       vac_tv_type = 'scenarios',
                                       vac_cov_ts = flu_vac_by_age,
                                       summary_type = 'diff',
                                       sim_startday = as.numeric(as.Date('2025-03-01') - as.Date('2020-01-01')) + 1,
                                       baseline_params = baseline_params,
                                       scale_us = T,
                                       scaled_metrics = primary_metrics)
  
  treat_vac_state_df$cumu_cases[i] = summary_traj['cumu_cases']
  treat_vac_state_df$cumu_hosp[i] = summary_traj['cumu_hosp']
  treat_vac_state_df$cumu_death[i] = summary_traj['cumu_death']
  
}

save(treat_vac_state_df,file = 'R/Rdata/baseline_by_state_548.Rdata')





