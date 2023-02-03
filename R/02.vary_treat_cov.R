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
library(dplyr)
library(ggplot2)
library(metR)

### load functions
source("R/cpp-v10/traj.from.params-v7-TT.R")
source("R/cpp-v10/main_traj.R")
source("R/functions/run_simulation_functions.R")
source("R/functions/generate_beta_func2.R")
source('R/functions/US_scale_up.R')

#### run simulation.R and baseline_params.func should be under the same dir ###
# load('R/cpp-v10/main/20230113/params/initial_params_548_80norebound.Rdata')
load('R/cpp-v10/main/20230113/params/initial_params_548.Rdata')

### load the tv-flu-vac-cov ###
load('data/us_flu_vac/relative_flu_vac_perc_by_age.Rdata')

############ Setting up baseline parameters ####################################

########## 3 transmission scenarios: half beta, beta, twice beta ###########
vac_cov = 0.49 # opt assumption

all_params = lapply(1:length(all_params_548),function(x) {
  
  ## add status quo vac cov ##
  vac_cover_generator(cover = vac_cov,
                      pop_by_age = get_pop_by_age('RI'),
                      vac_tv_type = 'scenarios',
                      vac_cov_ts = flu_vac_by_age,
                      baseline_params = all_params_548[[x]])
  
})
names(all_params) = c('opt','neu','pes')


############################# Simulation #######################################
### needed metrics ###
primary_metrics = c('cumu_cases','cumu_hosp','cumu_death')

##### different treat coverage run ###
# treat_cov = matrix(rep(c(0.137,seq(0.2,1,length.out = 9)),8),nrow = 10,ncol = 8)
treat_cov = matrix(rep(c(seq(0.5,0.6,0.01),seq(0.9,1,0.01)),8),
                   nrow = length(c(seq(0.5,0.6,0.01),seq(0.9,1,0.01))),
                   ncol = 8)

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
# save(all_treat_proj_df,file = 'R/cpp-v10/Rdata/20230113/treat_only_proj_548_80norebound.Rdata')

### save summaries ###
names(all_treat_summaries) = names(all_params)
all_treat_summaries_df = plyr::ldply(all_treat_summaries)
colnames(all_treat_summaries_df)[which(colnames(all_treat_summaries_df) == '.id')] = 'scene'

save(all_treat_summaries_df,file = 'R/cpp-v10/Rdata/20230113/treat_only_summaries_sim548_0.01granule.Rdata')

########################## Plotting ############################################
# load('R/cpp-v10/Rdata/20230113/treat_only_proj_365_80norebound.Rdata')
# treat_test = treat_cov[,1]
# 
# ## plotting new cases/hosp/deaths ##
# 
# all_treat_proj_df %>%
#   group_by(scene,treat) %>%
#   mutate(newcases = c(0,diff(cumu_cases,1)),
#          newhosp = c(0,diff(cumu_hosp,1)),
#          newdeaths = c(0,diff(cumu_death,1))) %>%
#   dplyr::select(scene,date,treat,newcases,newhosp,newdeaths) %>%
#   melt(id.vars = c('scene','date','treat')) %>%
#   mutate(variable = factor(variable,levels = c('newcases','newhosp','newdeaths'),
#                            labels = c('New Cases','New Hospitalizataions','New Deaths'))) %>%
#   ggplot() +
#   geom_line(aes(x = date, y = value,color = treat)) +
#   scale_color_brewer(palette = 'Spectral',name = 'Treatment Coverage(%)',
#                      labels = treat_test*100) +
#   scale_x_date(date_breaks = '1 year',date_labels = '%Y',expand = c(0.01,0)) +
#   facet_grid(variable ~ scene,scales = 'free') +
#   xlab('') +
#   ylab('') +
#   theme_bw() +
#   theme(text = element_text(size = 15),
#         axis.text = element_text(size = 13),
#         plot.margin = unit(c(0.1,0.1,1.2,0.1),'cm'),
#         legend.position = c(0.5,-0.07),
#         legend.direction = 'horizontal') + 
#   guides(colour = guide_legend(nrow = 1))
# ggsave(filename = paste0('R/cpp-v10/outputs/20230113/','incidence_treat_cov_seasAmp0.2_365_80norebound.jpg'),
#        width = 22,height = 10,
#        units = 'in')
# 
# 
# ####### plot the cumulative ######
# all_treat_proj_df %>%
#   melt(id.vars = c('scene','treat','date')) %>%
#   mutate(variable = factor(variable,levels = c('cumu_cases','cumu_hosp','cumu_death'),
#                            labels = c('Cases','Hospitalized','Deaths')),
#          scene = factor(scene,levels = c('Optimistic','Neutral','Pessimistic'),
#                         labels = c('Optimistic','Neutral','Pessimistic'))) %>%
#   ggplot() +
#   geom_line(aes(x = date,y = value,color = treat)) +
#   scale_color_brewer(palette = 'Spectral',name = 'Treatment Coverage (%)',
#                      labels = treat_test*100) +
#   scale_x_date(date_breaks = '1 year',date_labels = '%Y',
#                expand = c(0.01,0)) +
#   facet_grid(variable ~ scene,scales = 'free') +
#   xlab('') +
#   ylab('') +
#   theme_bw() +
#   theme(text = element_text(size = 15),
#         axis.text = element_text(size = 13),
#         plot.margin = unit(c(0.1,0.1,1.2,0.1),'cm'),
#         legend.position = c(0.5,-0.07),
#         legend.direction = 'horizontal') + 
#   guides(colour = guide_legend(nrow = 1))
# ggsave(filename = paste0('R/cpp-v10/outputs/20230113/','cumulative_treat_only_seasAmp0.2_365_80norebound.jpg'),
#        width = 22,height = 10,
#        units = 'in')


#### plot the data summaries during simulation under treatment ###
# load('R/cpp-v10/Rdata/20230113/treat_only_summaries_sim365_50TE.Rdata')
# 
# primary_summary_metrics = c('cumu_cases','cumu_hosp','cumu_death')
# all_treat_summaries_df$treat = as.character(all_treat_summaries_df$treat)
# all_treat_summaries_df$scene = factor(all_treat_summaries_df$scene, levels = c('opt','neu','pes'),
#                                labels = c('Optimistic','Neutral','Pessimistic'))
# ylabs = c('Total Cases','Total Hospitalizations','Total Deaths')
# palettes = c('Blues','Purples','YlOrBr')
# 
# 
# for(i in 1:length(primary_summary_metrics)) {
#   
#   ggplot(all_treat_summaries_df) +
#     geom_col(aes_string(x = 'scene',
#                         y = primary_summary_metrics[i],
#                         fill = 'treat'),
#              position = 'dodge') +
#     scale_fill_brewer(palette = palettes[i],
#                       name = 'Treatment\nCoverage(%)',
#                       labels = seq(0,1,length.out = 9)*100) +
#     xlab('') +
#     ylab(ylabs[i]) + 
#     theme_bw()
#   ggsave(filename = paste0('R/cpp-v10/outputs/20230113/sim_',primary_metrics[i],'_summary_treat_only_365_50TE.jpg'),
#          width = 8,height = 4,
#          units = 'in')
# }


