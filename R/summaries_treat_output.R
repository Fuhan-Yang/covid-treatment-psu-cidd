rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
### set up working dir to be github repo 
setwd('~/github/covid19-post-vaccination-burden')

load('R/cpp-v10/Rdata/20230113/treat_only_summaries_sim548_0.01granule.Rdata')
### annual numbers ###
years = 8

annual_all = 
  all_treat_summaries_df %>%
  mutate(cumu_cases = cumu_cases/(years * 1e6),
         cumu_hosp = cumu_hosp/(years * 1e3),
         cumu_death = cumu_death/(years * 1e3),
         treat = as.numeric(treat))

annual_all= annual_all[annual_all$treat > 0.13,]

annual_all$diff_cases = NA
annual_all$diff_hosp = NA
annual_all$diff_death = NA

annual_all$perc_cases = NA
annual_all$perc_hosp = NA
annual_all$perc_death = NA

new_annual_all = data.frame()

for(scene in unique(annual_all$scene)) {
  
  df = annual_all[annual_all$scene == scene, ]
  for(i in 1:nrow(df)) {
    
    df$diff_cases[i] = df$cumu_cases[df$treat == 0.137] - df$cumu_cases[i]
    df$perc_cases[i] = (df$cumu_cases[df$treat == 0.137] - df$cumu_cases[i])/df$cumu_cases[df$treat == 0.137] * 100
    
    df$diff_hosp[i] = df$cumu_hosp[df$treat == 0.137] - df$cumu_hosp[i]
    df$perc_hosp[i] = (df$cumu_hosp[df$treat == 0.137] - df$cumu_hosp[i])/df$cumu_hosp[df$treat == 0.137] * 100
    
    df$diff_death[i] = df$cumu_death[df$treat == 0.137] - df$cumu_death[i]
    df$perc_death[i] = (df$cumu_death[df$treat == 0.137] - df$cumu_death[i])/df$cumu_death[df$treat == 0.137] * 100
  }
   new_annual_all = rbind(new_annual_all,
                          df)
}

save(new_annual_all, file = 'R/cpp-v10/Rdata/20230113/treat_averted_summaries_548.Rdata')
load('R/cpp-v10/Rdata/20230113/treat_averted_summaries_548.Rdata')

new_annual_all %>%
  filter(treat == 1)%>%
  summarise(m_hosp = mean(perc_hosp),
            m_death = mean(perc_death),
            m_case = mean(perc_cases))

### reduction rate ###
red_r = 
  new_annual_all %>%
  group_by(scene) %>%
  mutate(r_case = c(0,diff(cumu_cases,1))*1e6,
         r_hosp = c(0,diff(cumu_hosp,1))*1000,
         r_death = c(0,diff(cumu_death,1))*1000)

# mean reduction rate on every 10% increase of treat cov #
red_r %>%
  group_by(scene) %>%
  filter(treat > 0.2) %>%
  summarise(case_m = mean(r_case),
            hosp_m = mean(r_hosp),
            death_m = mean(r_death)) 

#### plot the annual numbers #####
load('R/cpp-v10/Rdata/20230113/treat_averted_summaries_548.Rdata')

plot_annual_all = new_annual_all
plot_annual_all$treat = round(plot_annual_all$treat,3)

plot_annual_all = plot_annual_all[plot_annual_all$treat %in% c(0.137,0.2,0.4,0.6,0.8,1),]

primary_metrics = c('cumu_cases','cumu_hosp','cumu_death')
plot_annual_all$treat = as.character(plot_annual_all$treat)

plot_annual_all$scene = factor(plot_annual_all$scene, levels = c('opt','neu','pes'),
                               labels = c('Optimistic','Neutral','Pessimistic'))
ylabs = c('Annual Cases (millions)','Annual Hospitalized (thousands)','Annual Deaths (thousands)')



palettes = list(case = colorRampPalette(rev(brewer.pal(4,'Blues'))),
                hosp = colorRampPalette(rev(brewer.pal(4,'Purples'))),
                death = colorRampPalette(rev(brewer.pal(4,'YlOrBr'))))

for(i in 1:length(primary_metrics)) {
  
  p = ggplot() + 
    geom_col(aes_string(x = 'scene',
                        y = primary_metrics[i],
                        fill = 'treat'),
             position = 'dodge',
             data = plot_annual_all) +
    annotate('segment', x = 1.5,xend = 1.5, y = 0,yend = max(plot_annual_all[,primary_metrics[i]])*1.1,color = 'lightgrey')+ 
    annotate('segment', x = 2.5,xend = 2.5, y = 0,yend = max(plot_annual_all[,primary_metrics[i]])*1.1,color = 'lightgrey')+ 
    scale_x_discrete(expand = c(0.16,0.16),labels = NULL) + 
    scale_y_continuous(n.breaks = 10,expand = c(0,0),limits = c(0, max(plot_annual_all[,primary_metrics[i]])*1.1)) + 
    scale_fill_manual(values = palettes[[i]](length(unique(plot_annual_all$treat))),
                      name = 'Treatment\nCoverage(%)',
                      labels = as.numeric(unique(plot_annual_all$treat))*100) +
    # scale_fill_brewer('Blues',direction = -1) + 
    xlab('') +
    ylab(ylabs[i]) + 
    theme_bw() + 
    theme(panel.grid.major.x = element_blank(),
          text = element_text(size = 14),
          axis.text = element_text(size = 14))
  if(i == 3) {p = p + scale_x_discrete(expand = c(0.16,0.16))}
  ggsave(p, filename = paste0('R/cpp-v10/outputs/20230113/annual_',primary_metrics[i],'_summary_treat_only_548.jpg'),
         width = 6,height = 3,
         units = 'in')
}

######## plot averted burden ########
# plot_annual_all = new_annual_all
# plot_annual_all = plot_annual_all[-which(plot_annual_all$treat == 0.137),]
# 
# primary_summary_metrics = c('diff_cases','diff_hosp','diff_death')
# plot_annual_all$treat = as.character(plot_annual_all$treat)
# plot_annual_all$scene = factor(plot_annual_all$scene, levels = c('opt','neu','pes'),
#                                       labels = c('Optimistic','Neutral','Pessimistic'))
# ylabs = c('Annual Cases (millions)','Annual Hospitalizations (thousands)','Annual Deaths (thousands)')
# palettes = c('Blues','Purples','YlOrBr')
# 
# for(i in 1:length(primary_summary_metrics)) {
#   
#   ggplot(plot_annual_all) +
#     geom_col(aes_string(x = 'scene',
#                         y = primary_summary_metrics[i],
#                         fill = 'treat'),
#              position = 'dodge') +
#     scale_y_continuous(n.breaks = 10) + 
#     scale_fill_brewer(palette = palettes[i],
#                       name = 'Treatment\nCoverage(%)',
#                       labels = as.numeric(unique(plot_annual_all$treat))*100) +
#     xlab('') +
#     ylab(ylabs[i]) + 
#     theme_bw() 
#   ggsave(filename = paste0('R/cpp-v10/outputs/20230113/',primary_summary_metrics[i],'_summary_treat_only_548.jpg'),
#          width = 8,height = 4,
#          units = 'in')
# }


