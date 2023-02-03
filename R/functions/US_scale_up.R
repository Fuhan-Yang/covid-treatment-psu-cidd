library(dplyr)

###################### 2022-10-31 version: one ratio for one phase #############

get_ri_us_ratio = function(ri_fname = 'data/RI/RI-20221208.csv',
                           us_fname = 'data/US/data_table_for_weekly_case_trends__the_united_states.csv') {
  
  ri = read.csv(ri_fname)
  us = read.csv(us_fname,
                skip = 2)
  
  ### scale ri to us using the ratio of cases during four phases ###
  
  ### four phases: 
  ### before alpha: 2020-03-02 - 2021-03-29
  ### alpha: 2021-03-30 - 2021-06-22
  ### delta: 2021-06-23 - 2021-12-20
  ### omicron: 2021-12-21 - present 
  
  us_month = match(gsub('(\\D{3}).*','\\1',us$Date),month.abb)
  us_year = as.numeric(gsub('.*(\\d{4})$','\\1',us$Date))
  us_day = as.numeric(gsub('\\D{3}(.*)\\s\\d{4}$','\\1',us$Date))
  
  us_date = as.Date(paste0(us_year,'-',us_month,'-',us_day))
  us$date = us_date
  us = us[order(us$date),]
  
  ### Get the ratio from 2020-03-01 to 2022-11-30 ##
  us = us[us$date > as.Date('2020-02-29') & us$date < as.Date('2022-12-01'),]
  
  ri$date = as.Date(ri$date)
  ri = ri[ri$date > as.Date('2020-02-29') & ri$date < as.Date('2022-12-01'),]
  
  common_period = c(ri$date,us_date)[which(duplicated(c(ri$date,us_date)))]
  
  ### start after 2022-03-01 ###
  common_period = common_period[common_period > as.Date('2020-02-29')]
  ri = ri[ri$date %in% common_period,]
  ri$inc = c(0,diff(ri$cumulative_confirmed_cases,1))
  
  us = us[us$date %in% common_period,]
  
  ri_us_cases = data.frame(date = common_period, 
                           ri_inc = ri$inc,
                           us_inc = us$Weekly.Cases)
  
  ri_us_cases = ri_us_cases[order(ri_us_cases$date),]
  
  #### seperate the ri us data into four variant phases ####
  variant_date = c(as.Date('2020-03-01'),as.Date('2021-03-30'),as.Date('2021-06-23'),as.Date('2021-12-21'))
  
  ri_us_cases$var = findInterval(ri_us_cases$date,variant_date)
  
  ri_us_cases = 
    ri_us_cases %>%
    group_by(var) %>%
    mutate(ratio = sum(us_inc)/sum(ri_inc))
  
  ratios = data.frame(date = variant_date,
                      ratio = unique(ri_us_cases$ratio))
  
  return(ratios)
}


####################### 2022-10-22 version deprecated ##########################
# get_ri_us_ratio = function(ri_fname = 'data/RI/RI-partialformatted-20221007.csv',
#                        us_fname = 'data/US/data_table_for_daily_case_trends__the_united_states.csv') {
#   
#   ri = read.csv(ri_fname)
#   us = read.csv(us_fname,
#                 skip = 2)
#   
#   ### scale ri to us using ratio of averaged monthly incidence ###
#   ### get monthly incidence ###
#   ri$date = as.Date(ri$date)
#   ri_cases = data.frame(date = seq(min(ri$date),max(ri$date),'30 days'))
#   ri_cases$cumu = ri$cumulative_confirmed_cases[which(ri$date %in% ri_cases$date)]
#   
#   us_month = match(gsub('(\\D{3}).*','\\1',us$Date),month.abb)
#   us_year = as.numeric(gsub('.*(\\d{4})$','\\1',us$Date))
#   us_day = as.numeric(gsub('\\D{3}(.*)\\s\\d{4}$','\\1',us$Date))
#   
#   us_date = as.Date(paste0(us_year,'-',us_month,'-',us_day))
#   us$date = us_date
#   us = us[order(us$date),]
#   
#   us$cumu = cumsum(us$New.Cases)
#   
#   common_period = c(ri$date,us_date)[which(duplicated(c(ri$date,us_date)))]
#   common_month = seq(as.Date('2020-03-01'),as.Date('2022-10-01'),'month')
#   
#   ri_monthly_cases = data.frame(date = common_month,
#                                 cumu = ri$cumulative_confirmed_cases[ri$date %in% common_month])
#   
#   ri_monthly_cases$inc = c(0,diff(ri_monthly_cases$cumu,1))
#   
#   us_monthly_cases = data.frame(date = common_month,
#                                 cumu = us$cumu[which(us$date %in% common_month)])
#   
#   us_monthly_cases$inc = c(0,diff(us_monthly_cases$cumu,1))
#                               
#   us_monthly_cases$cumu = cumsum(us_monthly_cases$inc)
# 
#   #### Get the mean inc ratio for each month across 2020-2022 ###
#   # ratios = data.frame(date = common_month,
#   #                     ratio = us_monthly_cases$cumu/ri_monthly_cases$cumu)
#   
#   ## the ratio between inc varies month by month ##
#   ratios = data.frame(date = common_month,
#                       ratio = us_monthly_cases$inc/ri_monthly_cases$inc)
#   
#   ### delete 2020-03 for under-reporting ###
#   ratios = ratios[-which(ratios$date == as.Date('2020-03-01')),]
#   
#   ratios$month = format(ratios$date,'%m')
#   
#   m_ratios = 
#     ratios %>%
#     dplyr::group_by(month) %>%
#     dplyr::summarise(m = mean(ratio))
#   
#   return(m_ratios)
# }

#### Use the monthly ratio of cumulative cases between US and RI to scale up to the 
### rate of change in the US in terms of cases, hosps, and deaths
#'@param ri_proj, simulation from RI
#'@param ratio, ratio of the total cases between US and RI for four phases
#'@param metrics metrics that need to scale up to US, should be cumulative
#'@return scaled metrics of US

scale_to_us = function(ri_proj, 
                       ratios,
                       metrics = c('cumu_cases','cumu_hosp','cumu_death')) {
  
  ri_proj = data.frame(ri_proj)
  data_endday = as.Date('2020-02-29')
  
  ri_sim = ri_proj[ri_proj$date > data_endday,]
  
  ri_sim = ri_sim[,which(colnames(ri_sim) %in% c('date',metrics))]
  
  ri_sim_rate = data.frame(date = ri_sim$date)
  
  for(metric in metrics) {
    
    ri_sim_rate[,metric] = c(0,diff(ri_sim[,metric],1))
  }
  
  ri_sim_rate$ratio = ratios$ratio[findInterval(ri_sim_rate$date,ratios$date)]
  
  us_sim_rate = ri_sim_rate
  us_sim = data.frame(date = us_sim_rate$date)

  for(metric in metrics) {
    
   us_sim_rate[,metric] = ri_sim_rate[,metric] * ri_sim_rate$ratio
   us_sim[,metric] = cumsum(us_sim_rate[,metric])
  }
  
  return(us_sim)
  
}


