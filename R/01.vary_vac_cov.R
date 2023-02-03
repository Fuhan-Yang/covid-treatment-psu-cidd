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
library(RColorBrewer)
# library(brewer.pal)
library(metR)

### load functions
source("R/cpp-v10/traj.from.params-v7-TT.R")
source("R/cpp-v10/main_traj.R")
source("R/functions/run_simulation_functions.R")
source("R/functions/generate_beta_func2.R")
source('R/functions/US_scale_up.R')

#### run simulation.R and baseline_params.func should be under the same dir ###

### load params under three scenarios ###
# load('R/cpp-v10/main/20230113/params/initial_params_365.Rdata')
load('R/cpp-v10/main/20230113/params/initial_params_548_88VE.Rdata')
# load('R/cpp-v10/main/20230113/params/initial_params_730.Rdata')

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

all_params = lapply(1:length(all_params_548),function(x) {
  
  param_changer(change_params = treat_cover,
                baseline_params = all_params_548[[x]])
  
})
names(all_params) = c('opt','neu','pes')

############### Simulation #####################################################
### needed metrics ###
primary_metrics = c('cumu_cases','cumu_hosp','cumu_death')

##### different vac coverage run ###
# vac_test = c(0.49,seq(0.6,1,0.05))
# vac_test = c(seq(0.61,0.64,0.01),seq(0.81,0.89,0.01),seq(0.96,0.99,0.01))
# vac_test = seq(0.45,0.55,0.01)
vac_test = seq(0.66,0.7,0.01)
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

# save(all_vac_proj,file = 'R/cpp-v10/Rdata/20230113/vac_only_proj_365_0.01granule4.Rdata')

### save summaries ###
names(all_vac_summaries) = names(all_params)
all_vac_summaries = plyr::ldply(all_vac_summaries)

save(all_vac_summaries,file = 'R/cpp-v10/Rdata/20230113/vac_only_summaries_sim548_88VE_0.66-0.7.Rdata')

########################## Plotting ############################################

## plotting new cases/hosp/deaths ##
# load('R/cpp-v10/Rdata/20230113/vac_only_proj_548_50VE.Rdata')
# 
# all_new_proj =
#   all_vac_proj %>%
#   group_by(.id,vac) %>%
#   mutate(newcases = c(0,diff(cumu_cases,1)),
#          newhosp = c(0,diff(cumu_hosp,1)),
#          newdeaths = c(0,diff(cumu_death,1))) %>%
#   dplyr::select(.id,date,vac,newcases,newhosp,newdeaths)
# new_metrics = c('newcases','newhosp','newdeaths')
# 
# ylabs = c('New Cases','New Hospitalizataions','New Deaths')
# 
# molten_df = melt(all_new_proj,id.vars = c('.id','date','vac'))
# molten_df$variable = factor(molten_df$variable,levels = new_metrics,
#                             labels = c('New Cases','New Hospitalizataions','New Deaths'))
# 
# ### Add color for each class ###
# ggplot(molten_df) +
#   geom_line(aes(x = date, y = value,color = vac)) +
#   scale_color_brewer(palette = 'Spectral',name = 'Vaccination Coverage(%)',
#                      # limits = as.character(seq(0,1,0.1)),
#                      labels = vac_test*100) +
#   scale_x_date(date_breaks = '1 year',date_labels = '%Y',
#                expand = c(0,0)) +
#   xlab('') +
#   ylab('') +
#   facet_grid(variable ~ .id,scales = 'free') + 
#   theme_bw() +
#   theme(text = element_text(size = 15),
#         axis.text = element_text(size = 13),
#         plot.margin = unit(c(0.1,0.1,1.2,0.1),'cm'),
#         legend.position = c(0.5,-0.07),
#         legend.direction = 'horizontal') +
#   guides(colour = guide_legend(nrow = 1))
# ggsave(filename = paste0('R/cpp-v10/outputs/20230113/','incidence_vac_cov_seasAmp0.2_548_50VE.jpg'),
#        width = 22,height = 10,
#        units = 'in')
# 
# #### plot the vac only trajs ####
# ylabs = c('Total Cases','Total Hospitalizations','Total Deaths')
# palettes = c('Blues','Purples','YlOrBr')
# 
# molten_vac_proj = melt(all_vac_proj,id.vars = c('.id','date','vac'))
# 
# all_vac_proj %>%
#   melt(id.vars = c('.id','date','vac')) %>%
#   rename(scene = '.id') %>%
#   mutate(variable = factor(variable, levels = primary_metrics,
#                            labels = c('Total Cases','Total Hospitalizations','Total Deaths'))) %>%
#   ggplot() + 
#   geom_line(aes(x = date, y = value, color = vac)) + 
#   scale_x_date(date_breaks = '1 year',date_labels = '%Y') +
#   scale_color_brewer(palette = 'Spectral',
#                      name = 'Vaccination Coverage(%)',
#                      labels = vac_test*100) +
#   facet_grid(variable ~ scene,scales = 'free') +
#   xlab('') + 
#   ylab("") + 
#   theme_bw() + 
#   theme(text = element_text(size = 15),
#         axis.text = element_text(size = 13), 
#         plot.margin = unit(c(0.1,0.1,1.2,0.1),'cm'),
#         legend.position = c(0.5,-0.07),
#         legend.direction = 'horizontal') + 
#   guides(colour = guide_legend(nrow = 1))
# ggsave(filename = paste0('R/cpp-v10/outputs/20230113/cumulative_vac_proj_seasAmp0.2_548_50VE.jpg'),
#        width = 22,height = 10,
#        units = 'in')
  
#### plot the data summaries during simulation under vac only ###
# load('R/cpp-v10/Rdata/20230113/vac_only_summaries_sim548.Rdata')
# 
# primary_summary_metrics = c('cumu_cases','cumu_hosp','cumu_death')
# all_vac_summaries$vac = as.character(all_vac_summaries$vac)
# all_vac_summaries$.id = factor(all_vac_summaries$.id, levels = c('opt','neu','pes'),
#                           labels = c('Optimistic','Neutral','Pessimistic'))
# ylabs = c('Total Cases','Total Hospitalizations','Total Deaths')
# palettes = c('Blues','Purples','YlOrBr')
# 
# for(i in 1:length(primary_summary_metrics)) {
#   
#   ggplot(all_vac_summaries) +
#     geom_col(aes_string(x = '.id',
#                         y = primary_summary_metrics[i],
#                         fill = 'vac'),
#              position = 'dodge') +
#     scale_fill_brewer(palette = palettes[i],
#                       name = 'Vaccination\nCoverage(%)',
#                       labels = vac_test*100) +
#     xlab('') +
#     ylab(ylabs[i]) + 
#     theme_bw()
#   ggsave(filename = paste0('R/cpp-v10/outputs/20230113/sim_',primary_metrics[i],'_summary_vac_only_548.jpg'),
#          width = 8,height = 4,
#          units = 'in')
# }


