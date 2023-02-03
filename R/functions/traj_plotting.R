# setwd("~/github/covid19-post-vaccination-burden/cpp-v8")
# out = read.delim(file = 'ri.out')
# setwd("~/Dropbox/post_vaccination_covid_burden/for Fuhan")
# data = read.csv('RI_formatted_20210606-TT.csv')

plot_traj = function(out,data, plot_fname="traj_plot.jpg"){
  
  #### adjust the time frame of simulation ###
  # out = out[out[,1] %in% data$daynum,]
  
  #### add date, daynum = 1 is 2020-01-01 ###
  startdate = as.Date('2020-01-01')
  data$date = startdate + data$daynum - 1
  
  ### Save sim to plot in one data frame ###
  # out_to_plot = data.frame(date = data$date)
  out_to_plot = data.frame(date = seq(startdate + 49, startdate + out[nrow(out), 1] - 1, "1 day"))
  data_to_plot = data.frame(date = data$date)
  
  #### New Cases ###
  current_cases_cumu = rowSums(out[1:(nrow(out)-1),(290:298)])
  nextday_cases_cumu = rowSums(out[2:nrow(out),(290:298)])
  
  new_cases = nextday_cases_cumu - current_cases_cumu
  out_to_plot$new_cases = c(new_cases,NA)
  
  new_obs_cases = data$cumulative_confirmed_cases[2:nrow(data)] - data$cumulative_confirmed_cases[1:(nrow(data)-1)]
  data_to_plot$new_cases = c(new_obs_cases,NA)
  
  #### New Hospitalized ###
  current_hosp_cumu = rowSums(out[1:(nrow(out) -1),(299:307)])
  nextday_hosp_cumu = rowSums((out[2:nrow(out),(299:307)]))
  
  new_hosp = nextday_hosp_cumu - current_hosp_cumu
  out_to_plot$new_hosp = c(new_hosp,NA)
  
  new_obs_hosp = data$hospitalized_cumulative[2:nrow(data)] - data$hospitalized_cumulative[1:(nrow(data)-1)]
  data_to_plot$new_hosp = c(new_obs_hosp,NA)
  
  ### Current Hospitalized = HA +CA+CR+HR+V###
  HA = rowSums(out[,(137:172)])
  CA = rowSums((out[,(173:181)]))
  V = rowSums(out[,(182:235)])
  CR = rowSums(out[,(236:244)])
  HR = rowSums(out[,(245:253)])
  
  out_to_plot$current_hosp = rowSums(cbind(HA,CA,V,CR,HR))
  
  data_to_plot$current_hosp = data$hospitalized_currently
  
  ### New Death at home ###
  current_death_cumu = rowSums(out[1:(nrow(out)-1),(254:262)])
  nextday_death_cumu = rowSums(out[2:(nrow(out)),(254:262)])
  
  new_death = nextday_death_cumu - current_death_cumu
  
  out_to_plot$new_homedeath = c(new_death,NA)
  
  new_obs_death = data$cumulative_deaths[2:nrow(data)] - data$cumulative_deaths[1:(nrow(data) - 1)]
  data_to_plot$new_homedeath = c(new_obs_death,NA)
  
  ### New Hospital Death ###
  current_hd_cumu = rowSums(out[1:(nrow(out)-1),(263:271)])
  nextday_hd_cumu = rowSums(out[2:nrow(out),(263:271)])
  
  new_hosp_death = nextday_hd_cumu - current_hd_cumu
  out_to_plot$new_hospdeath = c(new_hosp_death,NA)

  new_obs_hd = data$Cumulative_hospital_deaths[2:nrow(data)] - data$Cumulative_hospital_deaths[1:(nrow(data) - 1)]
  data_to_plot$new_hospdeath = c(new_obs_hd,NA)
  
  #### New Hospital Discharge ###
  # current_disc_cumu = rowSums(out[1:(nrow(out)-1),(281:289)])
  # nextday_disc_cumu = rowSums(out[2:nrow(out),(281:289)])
  # 
  # new_hosp_disc = nextday_disc_cumu - current_disc_cumu
  # out_to_plot$new_hospdisc = c(new_hosp_disc,NA)
  # 
  # new_obs_disc = data$Cumulative_hospital_discharges[2:nrow(data)] - data$Cumulative_hospital_discharges[1:(nrow(data) - 1)]
  # data_to_plot$new_hospdisc = c(new_obs_disc,NA)
  
  #### Current ICU ###
  out_to_plot$current_ICU = CA + CR + V
  data_to_plot$current_ICU= data$InICU_currently
  
  ### Current Ventilation ###
  out_to_plot$current_vent = V
  data_to_plot$current_vent = data$OnVentilator_Currently
  
  
  ########## Plotting ########### 
  # jpeg(filename = 'traj_plot.jpg',width = 800,height = 1000)
  jpeg(filename = plot_fname, width = 800,height = 1000)
  par(mfrow = c(4,2),mar = c(2,3,2,2))
  for(i in 2:ncol(data_to_plot)) {
    
    plot(out_to_plot$date,out_to_plot[,i],type = 'l',
          col = 'grey',
         main = colnames(data_to_plot)[i],
         xlab = '',ylab = '',
         xaxt = 'n',
         yaxt = 'n',
         cex.main = 2)
    axis(1, at = seq(min(out_to_plot$date),max(out_to_plot$date),'3 months'),
         labels = format(seq(min(out_to_plot$date),max(out_to_plot$date),'3 months'),'%Y-%m'),
         cex.axis = 2)
    axis(2,cex.axis = 2)
    points(data_to_plot$date,data_to_plot[,i],pch = 20,cex = 0.3)
    
  }
  dev.off()
  
  
}
  
# plot_traj(out = out, data = data)
