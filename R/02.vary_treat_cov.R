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
load('R/params//initial_params_365.Rdata')
# load('R/params/initial_params_548.Rdata')
# load('R/params/initial_params_730.Rdata')

### load the tv-flu-vac-cov ###
load('data/us_flu_vac/relative_flu_vac_perc_by_age.Rdata')

############ Setting up baseline parameters ####################################

########## 3 transmission scenarios: half beta, beta, twice beta ###########
vac_cov = 0.49 # opt assumption

all_params = lapply(1:length(all_params_365),function(x) {
  
  ## add status quo vac cov ##
  vac_cover_generator(cover = vac_cov,
                      pop_by_age = get_pop_by_age('RI'),
                      vac_tv_type = 'scenarios',
                      vac_cov_ts = flu_vac_by_age,
                      baseline_params = all_params_365[[x]])
  
})
names(all_params) = c('opt','neu','pes')


############################# Simulation #######################################
### needed metrics ###
primary_metrics = c('cumu_cases','cumu_hosp','cumu_death')

##### different treat coverage run ###
treat_cov = matrix(rep(c(0.137,seq(0.2,1,length.out = 9)),8),nrow = 10,ncol = 8)

colnames(treat_cov) = paste0('tv-treat-cov-',seq(10,80,10),'_2')

all_treat_proj = list()
all_treat_summaries = list()

#### different treat under different beta scenarios ####
for(i in 1:length(all_params)) {
  
  treat_traj = run_multiple_params(change_params = treat_cov,
                                   baseline_params = all_params[[i]],
                                   output_type = 'primary',
                                   scale_us = T)
  names(treat_traj) = treat_cov[,1]
  treat_traj_df = plyr::ldply(treat_traj)
  colnames(treat_traj_df)[which(colnames(treat_traj_df) == '.id')] = 'treat'
  
  treat_summary = run_multiple_params(change_params = treat_cov,
                                      baseline_params = all_params[[i]],
                                      output_type = 'summary',
                                      summary_type = 'diff',
                                      sim_startday = as.numeric(as.Date('2025-03-01') - as.Date('2020-01-01')) + 1,
                                      scale_us = T)
  names(treat_summary) = treat_cov[,1]
  treat_summary_df = plyr::ldply(treat_summary)
  colnames(treat_summary_df)[which(colnames(treat_summary_df) == '.id')] = 'treat'
  
  
  all_treat_proj[[i]] = treat_traj_df
  all_treat_summaries[[i]] = treat_summary_df
}


################ Saving outputs ################################################
# names(all_treat_proj) = names(all_params)
# all_treat_proj_df = plyr::ldply(all_treat_proj)
# colnames(all_treat_proj_df)[which(colnames(all_treat_proj_df) == '.id')] = 'scene'
# 
# ### convert treat to string to use it as facet
# all_treat_proj_df$treat = as.character(all_treat_proj_df$treat) # coverage 
# all_treat_proj_df$scene = factor(all_treat_proj_df$scene, levels = c('opt','neu','pes'),
#                                  labels = c('Optimistic','Neutral','Pessimistic'))
# save(all_treat_proj_df,file = 'R/cpp-v10/Rdata/20230113/treat_only_proj_365.Rdata')

### save summaries ###
names(all_treat_summaries) = names(all_params)
all_treat_summaries_df = plyr::ldply(all_treat_summaries)
colnames(all_treat_summaries_df)[which(colnames(all_treat_summaries_df) == '.id')] = 'scene'

save(all_treat_summaries_df,file = 'R/Rdata/treat_only_summaries_sim365.Rdata')


