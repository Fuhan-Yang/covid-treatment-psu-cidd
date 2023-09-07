rm(list = ls())
setwd("~/github/covid-treatment-psu-cidd/")
library(dplyr)
flu_raw = read.csv('data/us_flu_vac/Influenza_Vaccination_Coverage_for_All_Ages__6__Months_.csv')

flu = flu_raw %>%
  filter(Vaccine == 'Seasonal Influenza',
         Geography.Type == 'HHS Regions/National',
         Dimension.Type == 'Age') 

us_flu = flu %>% 
  filter(Geography == 'United States') %>%
  mutate(date = NA,
         estimate = as.numeric(Estimate....))

#### Modify the date to align with inter-year flu season
# flu season is from year1-07/08 to year2-05 
for(i in 1:nrow(us_flu)) {
  
  ### define the year ###
  year1 = gsub('(\\d{4}).*','\\1',us_flu$Season.Survey.Year[i])
  year2 = paste0('20',gsub('.*(\\d{2})','\\1',us_flu$Season.Survey.Year[i]))
  
  ## flu season in year 1
  if(us_flu$Month[i] < 6) {
    
    us_flu$date[i] = as.Date(paste0(year2,'-',us_flu$Month[i],'-01'))
  } 
  
  ## flu season in year 2
  if(us_flu$Month[i] > 6) {
    
    us_flu$date[i] = as.Date(paste0(year1,'-',us_flu$Month[i],'-01'))
    
  }
  
}

us_flu$date = as.Date(us_flu$date,origin = '1970-01-01')

#### Sanity check on population by age ###
us_flu %>%
  group_by(Dimension) %>%
  summarise(n = n(), m = mean(Sample.Size))

### Select exclusive age groups #
us_flu = us_flu %>%
  filter(Dimension %in% c('6 Months - 4 Years','5-12 Years',
                          '13-17 Years','18-49 Years','50-64 Years',
                          '≥65 Years'),
         
         ### exclude 2009 H1N1 ###
         date >= as.Date('2010-06-01')) %>%
  select(date, Dimension,estimate, Sample.Size, Season.Survey.Year) %>%
  rename(age = Dimension) %>%
  
  ### Add cumulative vaccinated population ###
  mutate(vac_count = round(estimate*Sample.Size))


library(ggplot2)  

ggplot(us_flu) +
  geom_line(aes(x = date,
                y = estimate,
                color = age)) + 
  scale_x_date(date_breaks = '1 year',date_labels = '%Y-%m')

ggplot(us_flu) +
  geom_line(aes(x = date,
                y = vac_count,
                color = age)) + 
  scale_x_date(date_breaks = '1 year',date_labels = '%Y-%m')


### add 0 to June ###
add0 = function(date, y) {
  
  recent_year = format(date[length(date)], '%Y')
  # df = data.frame(date, y)
  
  if(recent_year %in% c('2011','2012')) {
    
    last_year = as.character(as.numeric(recent_year) -1)
    df = data.frame(date = c(as.Date(paste0(last_year,'-07-01')),date),
                    y = c(0,y))
    
  } else {
    
    df = data.frame(date = date, 
                    y = y)
  }
  added = data.frame(date = as.Date(paste0(recent_year,'-06-01')),
                     y = 0)

  rbind(df, added)
  
}

#### linearly interpolate monthly data to weekly data ####
df_interpolate = function(date,y) {
  
  #weeks = seq(date[1],date[length(date)],'7 days')
  weeks = seq(date[1],as.Date('2021-01-01'),'7 days')
  pts = approx(x = date,y = y,xout = weeks)
  
  df = data.frame(weeks, 
                  y = pts$y)
  return(df)
}

### create whole year flu vac cov ###
vac_cov_generator = function(date, y) {
  
  added0_df = add0(date,y)
  df_interpolate(date = added0_df$date,
                 y = added0_df$y)
  
}


### Get the general yearly flu vac by taking the mean vac across all the flu season ###
## flu season is from August to July next year 
us_flu_by_age = split(us_flu,us_flu$age)

mean_vac_by_age = lapply(1:length(us_flu_by_age),function(x){
  
  ### split the vac by each season ###
  vac_by_season = split(us_flu_by_age[[x]],us_flu_by_age[[x]]$Season.Survey.Year)
  
  ### In each season, get the vac cov by week from the cumu vac ###
  vac_by_season = lapply(1:length(vac_by_season),function(i) {
    
    ### get the cov and count difference from cumu data ###
    vac_by_season[[i]]$vac_diff = diff(c(0,vac_by_season[[i]]$estimate))
    vac_by_season[[i]]$count_diff = diff(c(0,vac_by_season[[i]]$vac_count))
    
    ### add 0 for the unreporting months ###
    cov_df = add0(date = vac_by_season[[i]]$date,
                  y = vac_by_season[[i]]$vac_diff)
    count_df = add0(date = vac_by_season[[i]]$date,
                      y = vac_by_season[[i]]$count_diff)
    
    vac_cov_count_df = cbind(cov_df,count_df$y)
    
    colnames(vac_cov_count_df) = c('weeks','cov','count')
    return(vac_cov_count_df)
    
  })
  
  vac_by_season = do.call(rbind, vac_by_season)
  vac_by_season$month_day = format(vac_by_season$weeks,'%m-%d')
  
  ### Get the weekly mean vac across seasons ###
  mean_vac_cov = vac_by_season %>%
    group_by(month_day) %>%
    summarise(m_cov = mean(cov),
              m_count = round(mean(count)))
  
  ### add fake year ###
  mean_vac_cov$date = as.Date(paste0('2020-',mean_vac_cov$month_day))
  
  ### add fake next year to get whole year interpolated ##
  mean_vac_cov = rbind(mean_vac_cov,
                       data.frame(month_day = '01-01',
                                  m_cov = mean_vac_cov$m_cov[mean_vac_cov$date == as.Date('2020-01-01')],
                                  m_count = mean_vac_cov$m_count[mean_vac_cov$date == as.Date('2020-01-01')],
                                  date = as.Date('2021-01-01')))
                                  
  
  ### linearly interpolate weekly data ###
  weekly_cov = df_interpolate(date = mean_vac_cov$date,
                            y = mean_vac_cov$m_cov)
  
  weekly_count = df_interpolate(date = mean_vac_cov$date,
                                y = mean_vac_cov$m_count)
  weekly_vac = cbind(weekly_cov,weekly_count$y)
  
  colnames(weekly_vac) = c('weeks','cov','count')
  weekly_vac$age = names(us_flu_by_age)[[x]]
  
  return(weekly_vac)
})

names(mean_vac_by_age) = names(us_flu_by_age)

#### plot the mean vac by age 
mean_vac_df = do.call(rbind, mean_vac_by_age)

ggplot(mean_vac_df) +
  geom_line(aes(x = weeks, y = cov,color = age))

save(mean_vac_df, 
     file = 'data/us_flu_vac/us_mean_interpolated_flu_vac_across2010_2021_by_age.Rdata')




