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

### assume the same status quo treat cov in the future ###
treat_cover = rep(0.137,8)
names(treat_cover) = c("tv-treat-cov-10_2",
                       "tv-treat-cov-20_2",
                       "tv-treat-cov-30_2",
                       "tv-treat-cov-40_2",
                       "tv-treat-cov-50_2",
                       "tv-treat-cov-60_2",
                       "tv-treat-cov-70_2",
                       "tv-treat-cov-80_2")

all_params = lapply(1:length(all_params_365),function(x) {
  
  param_changer(change_params = treat_cover,
                baseline_params = all_params_365[[x]])
  
})
names(all_params) = c('opt','neu','pes')

############### Simulation #####################################################
### needed metrics ###
primary_metrics = c('cumu_cases','cumu_hosp','cumu_death')

##### different vac coverage run ###
vac_test = c(0.49,seq(0.6,1,0.05))
all_vac_proj = list()
all_vac_summaries = list()

#### different vac under different beta scenarios ####
for(i in 1:length(all_params)) {
  
  primary_trajs = list()
  summary_trajs = data.frame()
  for(j in 1:length(vac_test)) {
    
    ### incorporate the yearly flu vac pattern ###
    vac_test_params = vac_cover_generator(cover = vac_test[j], 
                                          pop_by_age = get_pop_by_age(loc = 'RI'),
                                          vac_tv_type = 'scenarios',
                                          vac_cov_ts = flu_vac_by_age,
                                          baseline_params = all_params[[i]])
    
    vac_test_traj = traj.from.params(beta = vac_test_params$beta,
                                     params = vac_test_params$params,
                                     const.params = vac_test_params$const.params,
                                     non.odesim.params=vac_test_params$non.odesim.params,
                                     introday = vac_test_params$introday,
                                     tf=vac_test_params$tf,
                                     odepath=vac_test_params$odepath,
                                     loc=vac_test_params$loc,
                                     symp = vac_test_params$symp,
                                     mean.report.rate = vac_test_params$mean.report.rate)
    
    ### trajectories ###
    primary_traj = get_ode_main_traj(vac_test_traj,
                                     scale_up = T,
                                     scaled_metrics = c('cumu_cases',
                                                        'cumu_hosp',
                                                        'cumu_death'))
    ### totals of trajectories ###
    summary_traj = traj_summary(primary_traj,
                                sim_startday = as.numeric(as.Date('2025-03-01') - as.Date('2020-01-01')) + 1,
                                options = 'diff')
    
    
    primary_traj$vac = vac_test[j]
    primary_trajs[[j]] = primary_traj
    
    
    summary_traj = c(summary_traj,'vac' = vac_test[j])
    summary_trajs = rbind(summary_trajs,summary_traj)
    
  }
  primary_trajs = do.call(rbind, primary_trajs)
  all_vac_proj[[i]] = primary_trajs
  
  colnames(summary_trajs) = names(summary_traj)
  all_vac_summaries[[i]] = summary_trajs
}

################ Saving outputs ################################################
names(all_vac_proj) = names(all_params)
all_vac_proj = plyr::ldply(all_vac_proj)

### convert treat to string to use it as facet
all_vac_proj$vac = as.character(all_vac_proj$vac) # coverage 
all_vac_proj$.id = factor(all_vac_proj$.id, levels = c('opt','neu','pes'),
                          labels = c('Optimistic','Neutral','Pessimistic'))

# save(all_vac_proj,file = 'R/Rdata/vac_only_proj_365.Rdata')

### save summaries ###
names(all_vac_summaries) = names(all_params)
all_vac_summaries = plyr::ldply(all_vac_summaries)

save(all_vac_summaries,file = 'R/Rdata/vac_only_summaries_sim365.Rdata')



