rm(list = ls())
setwd("~/github/covid19-post-vaccination-burden")
library(dplyr)
library(data.table)
library(ggplot2)

### Get the US state names in R ###
treat_data = fread('data/US/cumu_therapeutics_HHS_20221211.csv')
data(state)

state_df = data.frame(state.abb = state.abb,
                      state.name = state.name)

treat_by_state = treat_data[Jurisdiction %in% state.name,]

##### Get the paxlovid dosages ######
pax_by_state = 
  treat_by_state[, `:=`(pax_delivered = gsub(',','',`Paxlovid  Delivered`),
                        pax_admin = gsub(',','',`Paxlovid Administered`))] %>%
  rename(state = Jurisdiction) %>%
  mutate(pax_delivered = as.numeric(pax_delivered),
         pax_admin = as.numeric(pax_admin)) %>%
  select(state, pax_delivered, pax_admin)

pax_by_state$state.abb = NA 
for(i in 1:nrow(pax_by_state)) {

  pax_by_state$state.abb[i] = state_df$state.abb[which(state_df$state.name == pax_by_state$state[i])]  
  
}

##### Get the cases when paxlovid is available ####
case_by_state = fread('data/US/Weekly_United_States_COVID-19_Cases_and_Deaths_by_State.csv')

case_by_state %>%
  ggplot() + 
  geom_path(aes(x = date_updated,y = tot_cases,color = state))

### The total cases after paxlovid is available: Jan 1 to Dec 11 2022
# Choose the closest date: Dec 30 2021 to Dec 8 2022
case_by_state = case_by_state[state %in% state.abb,]

diff_case_by_state = 
  case_by_state %>%
  filter(date_updated == as.Date('2022-12-08')) %>%
  group_by(state) %>%
  pull(tot_cases) - 
  
  case_by_state %>%
  group_by(state)%>%
  filter(date_updated == as.Date('2021-12-30')) %>%
  pull(tot_cases)

treat_cov_by_state_df = data.frame(state = unique(case_by_state$state),
                                   cases = diff_case_by_state)

treat_cov_by_state_df$pax = NA
for(i in 1:nrow(treat_cov_by_state_df)) {

  treat_cov_by_state_df$pax[i] = pax_by_state$pax_admin[which(pax_by_state$state.abb == treat_cov_by_state_df$state[i])]  
  
}

### get the baseline treatment coverage for each state ###
treat_cov_by_state_df = 
  treat_cov_by_state_df %>%
  mutate(cov = round(pax/cases,3)*100)

treat_cov_by_state_df = treat_cov_by_state_df[order(treat_cov_by_state_df$cov,decreasing = T),]
save(treat_cov_by_state_df,file = 'R/cpp-v10/Rdata/treat_cov_baseline_by_state.Rdata')

load('R/cpp-v10/Rdata/treat_cov_baseline_by_state.Rdata')
ggplot(treat_cov_by_state_df) + 
  geom_point(aes(x = reorder(state,-cov), y = cov),size = 3) + 
  scale_y_continuous(n.breaks = 10) + 
  xlab('State') + 
  ylab('Treatment Coverage (%)') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = 'figs/supp_fig_treat_cov_by_state.jpg',width = 8,height = 4,
       units = 'in')


### get the national treatment coverage ###
us = fread('data/US/data_table_for_weekly_case_trends__the_united_states.csv')

us_month = match(gsub('(\\D{3}).*','\\1',us$Date),month.abb)
us_year = as.numeric(gsub('.*(\\d{4})$','\\1',us$Date))
us_day = as.numeric(gsub('\\D{3}(.*)\\s\\d{4}$','\\1',us$Date))

us_date = as.Date(paste0(us_year,'-',us_month,'-',us_day))
us$date = us_date
us = us[order(us$date),]  

us = 
  us %>%
  rename(new_case = `Weekly Cases`) %>%
  mutate(cumu_cases = cumsum(new_case))

plot(us$date,us$cumu_cases,type = 'l')

## national coverage ##
total_pax = 
  treat_data %>%
  filter(Jurisdiction == 'Grand Total') %>%
  mutate(pax_admin = as.numeric(gsub(',','',`Paxlovid Administered`))) %>%
  pull(pax_admin) 

total_cases = us %>%
  filter(date == as.Date('2022-12-07')) %>%
  pull(cumu_cases) - 
  us %>%
  filter(date == as.Date('2021-12-29')) %>%
  pull(cumu_cases)

# 0.137
total_pax/total_cases

# national admined paxlovid > sum of each state 
total_pax
sum(pax_by_state$pax_admin)

