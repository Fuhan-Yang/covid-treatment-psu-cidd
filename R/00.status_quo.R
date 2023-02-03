## Under cpp-v10 (with V class), the status quo condition during 2022-03-01 to 2023-03-01 ##

###############################################################################
### 1. Preliminaries
###############################################################################
rm(list = ls())
### set up working dir to be github repo 
setwd('~/github/covid19-post-vaccination-burden')

### necessary packages
## library for spline-based expansions
library(fda)
## package for iSplines:
library(splines2)
## package for plotting
library(reshape2)
library(ggplot2)
library(metR)

### load functions
source("R/cpp-v10/traj.from.params-v7-TT.R")
source("R/cpp-v10/main_traj.R")
source("R/functions/run_simulation_functions.R")
source("R/functions/generate_beta_func2.R")
source('R/functions/US_scale_up.R')


#### run simulation.R and baseline_params.func should be under the same dir ###
# when dur = 365, add mean_beta argument #
source("R/cpp-v10/main/20230113/calibrate/baseline_parameters_RI365.R")
# source("R/cpp-v10/main/20230113/calibrate/baseline_parameters_RI548.R")
# source("R/cpp-v10/main/20230113/calibrate/baseline_parameters_RI730.R")

### load the tv-flu-vac-cov ###
load('data/us_flu_vac/relative_flu_vac_perc_by_age.Rdata')

odepath = 'cpp-v10-v-class/'

#### status quo parameters ###
# baseline beta is during entire Omicron: 2021/12/21 - 2022/11/30
base_beta_start = 851
base_period = 90

params <- generate_baseline_params_RI(odepath = odepath,
                                      treat_beginday = 731, # treatment starts on Jan 1 2022
                                      treat_endday = 1156,
                                      beta_base_period_start = base_beta_start,# baseline beta is during entire Omicron: 2021/12/21 - 2022/11/30
                                      beta_base_period = base_period,
                                      mean_beta = 0.4283 * 1.2,
                                      sim_time = 91, # simulate until 03/01/2023
                                      beta_transition_period = 30,
                                      beta_multiplier = 1,
                                      seas_amp = 0, # no simulated seasonal increase at current season
                                      seas_dur = 120, # High season lasts four months
                                      seas_start = 210) # High season: Oct - Feb

vac_cover = 0.49 # conservative assumption of national coverage: 36%; optimistic: 49%

#### Add current vaccine coverage: 36% from Mar'22 to Mar'23
params_vac = vac_cover_generator(cover = vac_cover,
                                      pop_by_age = get_pop_by_age(loc = 'RI'),
                                      vac_tv_type = 'scenarios',
                                      vac_cov_ts = flu_vac_by_age,
                                      baseline_params = params)

#### Add RI current treatment coverage: 21.9% from Jan 1 to Dec 2022 ###
## Add US current treatment coveraeg: 13.7%
treat_cover = rep(0.137,8)
names(treat_cover) = c("tv-treat-cov-10_1",
                       "tv-treat-cov-20_1",
                       "tv-treat-cov-30_1",
                       "tv-treat-cov-40_1",
                       "tv-treat-cov-50_1",
                       "tv-treat-cov-60_1",
                       "tv-treat-cov-70_1",
                       "tv-treat-cov-80_1")

current_params = param_changer(change_params = treat_cover,
                               baseline_params = params_vac)
current_params_548 = current_params

# save(current_params_548,file = 'R/cpp-v10/main/20230113/params/status_quo_params548.Rdata')

status_quo_traj = traj.from.params(beta = current_params$beta,
                                   params = current_params$params,
                                   const.params = current_params$const.params,
                                   non.odesim.params=current_params$non.odesim.params,
                                   introday = current_params$introday,
                                   tf=current_params$tf,
                                   odepath=current_params$odepath,
                                   loc=current_params$loc,
                                   symp = current_params$symp,
                                   mean.report.rate = current_params$mean.report.rate)

primary_traj = get_ode_main_traj(status_quo_traj,
                                 scale_up = T,
                                 scaled_metrics = c('cumu_cases',
                                                    'cumu_hosp',
                                                    'cumu_death'))

startdate = as.Date('2022-03-01')
enddate = as.Date('2023-02-28')

### Deaths in the current season ###
primary_traj$cumu_death[primary_traj$date == enddate] -
  primary_traj$cumu_death[primary_traj$date == startdate]

#### current season deaths: 95% VE; 49% vac cov; 13.7% treat cov ###
# 365 day: 170,763 (beta =  0.4283 * 1.2)
# 548 day: 165,921 (beta: 0.808) 170,662 (50% VE)
# 730 day: 168,198 (beta: 1.273)


###### generate three transmission scenarios #######
########## 3 transmission scenarios: half beta, beta, twice beta ###########

### baseline period of beta: May 1 - July 1 2022
base_beta_start = 851
base_period = 60

opt_beta = generate_beta2(start_beta = current_params_548$beta,
                          base_period_start = base_beta_start, 
                          base_period = base_period,
                          extend_period = 3654,
                          multiplier = 0.5,
                          amp = 0.2,
                          seas_dur = 123, 
                          seas_start = 210)

neutral_beta = generate_beta2(start_beta = current_params_548$beta,
                              base_period_start = base_beta_start, 
                              base_period = base_period,
                              extend_period = 3654,
                              multiplier = 1,
                              amp = 0.2,
                              seas_dur = 123, 
                              seas_start = 210)

pes_beta = generate_beta2(start_beta = current_params_548$beta,
                          base_period_start = base_beta_start, 
                          base_period = base_period,
                          extend_period = 3654,
                          multiplier = 2,
                          amp = 0.2,
                          seas_dur = 123, 
                          seas_start = 210)

all_betas = list(opt_beta,neutral_beta,pes_beta)

#### beta is calibrated for 2022-2023, the beta onwards is the average beta from this season ###

time = seq(as.Date('2020-03-01'),as.Date('2033-02-28'),'1 day')
beta_df = data.frame(time = time,
                     beta = neutral_beta)

### avg beta in 2022-2023 548-day: 0.808 ###
mean(beta_df$beta[beta_df$time %in% seq(as.Date('2022-03-01'),as.Date('2023-02-28'),'1 day')])

### avg beta in the future seasons ###
mean(beta_df$beta[beta_df$time %in% seq(as.Date('2023-03-01'),as.Date('2024-02-28'),'1 day')])

### future beta should be scaled up to the same level as in 2022-2023
ratio = mean(beta_df$beta[beta_df$time %in% seq(as.Date('2022-03-01'),as.Date('2023-02-28'),'1 day')])/mean(beta_df$beta[beta_df$time %in% seq(as.Date('2023-03-01'),as.Date('2024-02-28'),'1 day')])

### the index to scale up beta
ind = which(time == as.Date('2023-03-01'))

all_betas = lapply(1:length(all_betas),function(x) {
  
  ### scale future beta to the level of 2022-2023 ###
  for(i in 1:length(all_betas[[x]])) {all_betas[[x]][i] = ifelse(i < ind,
                                                                 all_betas[[x]][i],
                                                                 all_betas[[x]][i] * ratio)}
  return(all_betas[[x]])
  
})

#### generate initial params for each scenario, need to append future cov ###
all_params_548 = lapply(1:length(all_betas),function(x) {
  
  current_params_548$beta = all_betas[[x]]
  current_params_548$tf = length(all_betas[[x]]) + 60
  
  return(current_params_548)
  
})
names(all_params_548) = c('opt','neu','pes')

all_params_365 = all_params_548
save(all_params_365,file = 'R/cpp-v10/main/20230113/params/initial_params_365_80norebound.Rdata')

#### Projections under current strategies ####
neu_params = all_params_548$neu

future_treat_cov = rep(0.137,8)
names(future_treat_cov) = c("tv-treat-cov-10_2",
                            "tv-treat-cov-20_2",
                            "tv-treat-cov-30_2",
                            "tv-treat-cov-40_2",
                            "tv-treat-cov-50_2",
                            "tv-treat-cov-60_2",
                            "tv-treat-cov-70_2",
                            "tv-treat-cov-80_2")

neu_params = param_changer(change_params = future_treat_cov,
                           baseline_params = neu_params)

vac_cover = 0.49
neu_params = vac_cover_generator(cover = vac_cover,
                                 pop_by_age = get_pop_by_age('RI'),
                                 vac_tv_type = 'scenarios',
                                 vac_cov_ts = flu_vac_by_age,
                                 baseline_params = neu_params)

all_traj = traj.from.params(beta = neu_params$beta,
                                 params = neu_params$params,
                                 const.params = neu_params$const.params,
                                 non.odesim.params=neu_params$non.odesim.params,
                                 introday = neu_params$introday,
                                 tf=neu_params$tf,
                                 odepath=neu_params$odepath,
                                 loc=neu_params$loc,
                                 symp = neu_params$symp,
                                 mean.report.rate = neu_params$mean.report.rate)

### trajectories of RI ###
# primary_traj_ri = get_ode_main_traj(all_traj,
#                                  scale_up = F,
#                                  scaled_metrics = c('cumu_cases',
#                                                     'cumu_hosp',
#                                                     'cumu_death'))
# 
# sus_traj = data.frame(date = as.Date('2020-01-01') + all_traj[,1] - 1,
#                       sus = rowSums(all_traj[,c(2:10)]))
# 
# rec_traj = data.frame(date = as.Date('2020-01-01') + all_traj[,1] - 1,
#                       rec = rowSums(all_traj[,c(272:280)]))
# 
# rec_hosp_traj = data.frame(date = as.Date('2020-01-01') + all_traj[,1] - 1,
#                            rec = rowSums(all_traj[,c(281:289)]))
# 
# vac_traj = data.frame(date = as.Date('2020-01-01') + all_traj[,1] - 1,
#                       rec = rowSums(all_traj[,c(290:298)]))

### trajectories of US ###
primary_traj = get_ode_main_traj(all_traj,
                                 scale_up = T,
                                 scaled_metrics = c('cumu_cases',
                                                    'cumu_hosp',
                                                    'cumu_death'))
### totals of trajectories ###
summary_traj = traj_summary(primary_traj,
                            sim_startday = as.numeric(as.Date('2025-03-01') - as.Date('2020-01-01')) + 1,
                            options = 'diff')

summary_traj/8



