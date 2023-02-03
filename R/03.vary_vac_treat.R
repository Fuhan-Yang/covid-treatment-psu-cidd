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

### load the tv-flu-vac-cov ###
load('data/us_flu_vac/relative_flu_vac_perc_by_age.Rdata')

############ load initial parameters: params fixed before Mar 1 2023 ###########
load('R/cpp-v10/main/20230113/params/initial_params_365.Rdata')
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
# save(all_param_combs_list,file = 'R/cpp-v10/Rdata/20230113/vac_treat_summaries_sim365.Rdata')
save(all_param_combs_list,file = 'R/cpp-v10/Rdata/20230113/vac_treat_summaries_sim365_20250301.Rdata')

load('R/cpp-v10/Rdata/20230113/vac_treat_summaries_sim365_20250301.Rdata')
primary_metrics = c('cumu_cases','cumu_hosp','cumu_death')

### Get the total cases/hosp/deaths per year in thousands ###

## 2025 - 2033 ##
years = 8


all_param_combs_list = lapply(1:length(all_param_combs_list),function(x) {
  
  all_param_combs_list[[x]] %>%
    filter(treat != 0.137,vac != 0.49) %>% ## remove 0.49 and 0.137 for better plotting ##
    mutate(cumu_cases = round(cumu_cases/(years*1e6),1), # in millions
           cumu_hosp = round(cumu_hosp/(years*1e3),1),   # in thousands
           cumu_death = round(cumu_death/(years*1e3),1)) %>%  # in thousands
    dplyr::select(treat,vac,cumu_cases,cumu_hosp,cumu_death)
})
names(all_param_combs_list) = c('opt','neu','pes')

metric_list = vector('list',length = 3)
names(metric_list) = primary_metrics


legend_names = c('Cases','Hospitalized','Deaths')
units = c('millions','thousands','thousands')
palettes = c('Blues','Purples','YlOrBr')

for(i in 1:length(all_param_combs_list)) {
  
  for(metric in primary_metrics) {
    
    ggplot(all_param_combs_list[[i]]) + 
      geom_raster(aes_string(x = 'treat', y = 'vac', fill = metric)) + 
      geom_text(aes_string(x = 'treat', y = 'vac', label = metric)) + 
      scale_fill_distiller(palette = palettes[which(metric == primary_metrics)],direction = 1,
                           name = paste0('Annual ',
                                         legend_names[which(metric == primary_metrics)],
                                         ' (',
                                         units[which(metric == primary_metrics)],
                                         ')'),
                           n.breaks = 6) + 
      scale_x_continuous(breaks = seq(0,1,0.05),expand = c(0,0),labels = paste0(seq(0,1,0.05)*100)) + 
      scale_y_continuous(breaks = seq(0,1,0.05),expand = c(0,0),labels = paste0(seq(0,1,0.05)*100)) + 
      xlab('Treatment Coverage (%)') +
      ylab('Vaccination Coverage (%)') +
      guides(fill = guide_colorbar(direction = 'horizontal',
                                   barwidth = 10,barheight = 0.5)) + 
      theme_bw() + 
      theme(text = element_text(size = 13),
            plot.margin = unit(c(0.1,0.1,1.2,0.1),'cm'),
            legend.position = c(0.5,-0.13))
    ggsave(filename = paste0('R/cpp-v10/outputs/20230113/heatmap_95VE_365_202503_',
                             gsub('cumu_(.*)','\\1',metric),'_',
                             names(all_param_combs_list[i]),'.jpeg'),
           width = 10, height = 6,units = 'in')
    # ggsave(filename = paste0('figs/contour/Fig4.contour_', 
    #                          gsub('cumu_(.*)','\\1',metric),
    #                          '_varying_vac_treat_',
    #                          names(all_param_combs_list)[i],
    #                          '_sim365.jpeg'),
    #        width = 10, height = 6,units = 'in')

    
  }
  
}
  
