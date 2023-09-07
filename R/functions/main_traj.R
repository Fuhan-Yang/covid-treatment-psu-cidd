##### This version adapts to cpp-v10 ##########
#' function to get the primary trajectories from traj.from.params()
#' @out is the direct output from traj.from.params()
#' @scale_up if TRUE,scale up from RI to US
#' @scaled_metrics only applies when scale_up is TRUE, defines the metrics that needs to scale up
#' @data_endday only applies when scale_up is TRUE, the end day of observed data
#' @return main trajectories: cumu_cases, cumu_hosp, cumu_death, current_hosp,
#' current_ICU, current_vent
get_ode_main_traj = function(out,
                             scale_up = F,
                             scaled_metrics = c('cumu_cases',
                                                'cumu_hosp',
                                                'cumu_death')){
  
  #### adjust the time frame of simulation ###
  ### start date is 2020-01-01 ###
  startdate = as.Date('2020-01-01')
  
  ### The simulation starts at day 50, add 49 in the real date ##
  ### round the daynum since it can be decimal ###
  ### Save sim to plot in one data frame ###
  out_to_plot = data.frame(date = seq(startdate + 49, startdate + round(out[nrow(out), 1]) - 1, "1 day"))
  
  ### Cumulative Cases ###
  out_to_plot$cumu_cases = rowSums(out[1:nrow(out),(299:307)])
  
  ## Cumulative Hospitalized ##
  out_to_plot$cumu_hosp = rowSums(out[1:nrow(out),(308:316)])
  
  ### Cumu Death at Hospital ###
  cumu_hospdeath = rowSums(out[1:nrow(out),(263:271)])
  ### Cumu Death at Home ##
  cumu_homedeath = rowSums(out[1:nrow(out),(254:262)])
  
  ### Total cumu death ###
  out_to_plot$cumu_death = rowSums(cbind(cumu_homedeath,cumu_hospdeath))
  

############################ Deprecated 2022/10/11 #############################  
  ### Current Hospitalized = HA +CA+CR+HR+V ###
  # HA = rowSums(out[,(137:172)])
  # CA = rowSums((out[,(173:181)]))
  # V = rowSums(out[,(182:235)])
  # CR = rowSums(out[,(236:244)])
  # HR = rowSums(out[,(245:253)])
  # 
  # out_to_plot$current_hosp = rowSums(cbind(HA,CA,V,CR,HR))
  # 
  # #### Current ICU ###
  # out_to_plot$current_ICU = CA + CR + V
  # 
  # ### Current Ventilation ###
  # out_to_plot$current_vent = V
#################### Deprecated 2022/10/11 #####################################
  
  if(scale_up) {
    us_sim = scale_to_us(ri_proj = out_to_plot,
                         ratios = get_ri_us_ratio(),
                         metrics = scaled_metrics)
    return(us_sim)
  } else {
    
    return(out_to_plot)
  }
}

#' plot the main trajectories returned by get_ode_main_traj
#' @data the observed data from state DOH
#' @out the output from get_ode_main_traj
#' @fname specify the location to save the plot
plot_main_traj = function(data, out, fname,save = T) {
  
  # out = get_ode_main_traj(out)
  
  #### add date, daynum = 1 is 2020-01-01 ###
  startdate = as.Date('2020-01-01')
  data$date = startdate + data$daynum - 1
  
  data_to_plot = data.frame(date = data$date)
  
  ### cumulative cases ###
  data_to_plot$cumu_cases = data$cumulative_confirmed_cases
  
  ### cumulative hosp ###
  data_to_plot$cumu_hosp = data$hospitalized_cumulative
  
  ### cumulative homedeath ###
  cumu_homedeath = data$cumulative_deaths
  
  ### cumulative hospdeath ###
  cumu_hospdeath = data$Cumulative_hospital_deaths
  
  ### cumulative total death ###
  data_to_plot$cumu_death = rowSums(cbind(cumu_homedeath,
                                          cumu_hospdeath))
  
  ### current hosp ###
  data_to_plot$current_hosp = data$hospitalized_currently
  
  ### current ICU ###
  data_to_plot$current_ICU = data$InICU_currently
  
  ### current Vent ###
  data_to_plot$current_vent = data$OnVentilator_Currently
  
  if(save == T) {
    jpeg(filename = paste0(fname,'.jpg'), width = 800,height = 1000)
    par(mfrow = c(3,2),mar = c(2,3,2,2))
    for(i in 2:ncol(data_to_plot)) {
      
      plot(out$date,out[,i],type = 'l',
           col = 'grey',
           main = colnames(data_to_plot)[i],
           xlab = '',ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 2)
      axis(1, at = seq(min(out$date),max(out$date),'3 months'),
           labels = format(seq(min(out$date),max(out$date),'3 months'),'%Y-%m'),
           cex.axis = 2)
      axis(2,cex.axis = 2)
      points(data_to_plot$date,data_to_plot[,i],pch = 20,cex = 0.3)
      
    }
    dev.off()
    
  } else {
    
    par(mfrow = c(3,2),mar = c(2,3,2,2))
    for(i in 2:ncol(data_to_plot)) {
      
      plot(out$date,out[,i],type = 'l',
           col = 'grey',
           main = colnames(data_to_plot)[i],
           xlab = '',ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 2)
      axis(1, at = seq(min(out$date),max(out$date),'3 months'),
           labels = format(seq(min(out$date),max(out$date),'3 months'),'%Y-%m'),
           cex.axis = 2)
      axis(2,cex.axis = 2)
      points(data_to_plot$date,data_to_plot[,i],pch = 20,cex = 0.3)
      
    }
  }
}


###### This function extracts the age-specific current vaccinated population ###
#'@param out output from traj.from.params
#'@scale_up logical, if scale up RI results to US 
#'@return data frame of age-stratified current vac data with date
get_vac_traj = function(out,scale_up = T) {
  
  #### adjust the time frame of simulation ###
  ### start date is 2020-01-01 ###
  startdate = as.Date('2020-01-01')
  
  ### The simulation starts at day 50, add 49 in the real date ##
  ### round the daynum since it can be decimal ###
  ### Save sim to plot in one data frame ###
  
  age_specific_out = data.frame(date = seq(startdate + 49, startdate + round(out[nrow(out), 1]) - 1, "1 day"),
                                age0 = NA,
                                age10 = NA,
                                age20 = NA,
                                age30 = NA,
                                age40 = NA,
                                age50 = NA,
                                age60 = NA,
                                age70 = NA,
                                age80 = NA)
  
  age_specific_out[,c(2:10)] = out[1:nrow(out),c(290:298)]
  
  if(scale_up) {
    
    us_sim = scale_to_us(ri_proj = age_specific_out,
                         ratios = get_ri_us_ratio(),
                         metrics = paste0('age',seq(0,80,10)))
    return(us_sim)
  } else {return(age_specific_out)}
  
}


### This function aims to extract the age-specific cumu cases, hosp, and deaths traj
### from traj.from.params ###

#'@param out output from traj.from.params 
#'@param metric choose the return metric from cases, hosp, and deaths
#'@return age-stratified cumulative metrics 

get_age_traj = function(out, 
                        metric,
                        scale_up = F) {
  
  #### adjust the time frame of simulation ###
  ### start date is 2020-01-01 ###
  startdate = as.Date('2020-01-01')
  
  ### The simulation starts at day 50, add 49 in the real date ##
  ### round the daynum since it can be decimal ###
  ### Save sim to plot in one data frame ###
  
  age_specific_out = data.frame(date = seq(startdate + 49, startdate + round(out[nrow(out), 1]) - 1, "1 day"),
                                age0 = NA,
                                age10 = NA,
                                age20 = NA,
                                age30 = NA,
                                age40 = NA,
                                age50 = NA,
                                age60 = NA,
                                age70 = NA,
                                age80 = NA)
  
  if(metric == 'cases') {age_specific_out[,c(2:10)] = out[1:nrow(out),(299:307)]}
  else if(metric == 'hosp') {age_specific_out[,c(2:10)] = out[1:nrow(out),(308:316)]}
  else if(metric == 'deaths') {
    
    ### Cumu Death at Hospital ###
    cumu_hospdeath = out[1:nrow(out),(263:271)]
    ### Cumu Death at Home ##
    cumu_homedeath = out[1:nrow(out),(254:262)]
    
    age_specific_out[,c(2:10)] = cumu_hospdeath + cumu_homedeath
  }
  else {stop('only cases,hosp, or deaths are available metrics,try another one.')}
  
  if(scale_up) {
    
    us_sim = scale_to_us(ri_proj = age_specific_out,
                         ratios = get_ri_us_ratio(),
                         metrics = paste0('age',seq(0,80,10)))
    return(us_sim)
  } else {return(age_specific_out)}
  
}


#' get the final output for cumulative metrics, get the max/mean/median current
#' metrics(deprecated at 09/15/2022)
#' @param out is the primary output from get_ode_main_traj
#' @param sim_startday the daynum that starts simulation, only needed when options == 'diff'
#' @param options defines summary type,
#' 'max': return the final/max of the cumulative metrics
#' 'diff': return the difference in the cumulative metrics after and before the estimated date
#' @return the final values for the cumulative metrics, 
#' and the max/mean/median/all for the current metrics
traj_summary = function(out,sim_startday,options) {
  
  metric_df = out[,-which(colnames(out) == 'date')]
  max_df = sapply(metric_df, max,na.rm = T)
  
  ## return the final/max values ##
  if(options == 'max') {
    
    return(max_df)
    
    ## return the difference between final values and the values before simulation ##
  } else if(options == 'diff') {
    
    ## default startday ##
    startday = as.Date('2020-01-01')
    sim_startdate = startday + sim_startday - 1
    
    ### get the difference between projection and current status ##
    est = out[out$date <= sim_startdate,]
    est_cumu_df = est[,-which(colnames(est) == 'date')]
    
    est_summ_df = sapply(est_cumu_df,max,na.rm = T)
    
    summ_df = max_df - est_summ_df
    
    return(summ_df)
  }
 
}

  
  
  ## current_temp is deprecated 2022/10/09 ## 
  ## for current metrics, options are max,mean,median
  # current_df = out[,grep('current.*',colnames(out))]
  # 
  # if(options == 'max') {
  #   current_temp = sapply(current_df,max,na.rm = T) 
  #   names(current_temp) = paste0(names(current_temp),'_max')
  # }
  # if(options == 'mean') {
  #   current_temp = sapply(current_df,mean,na.rm = T) 
  #   names(current_temp) = paste0(names(current_temp),'_mean')
  # }
  # if(options == 'median') {
  #   current_temp = sapply(current_df,median,na.rm = T) 
  #   names(current_temp) = paste0(names(current_temp),'_median')
  # }
  # if(options == 'all') {
  #   current_max = sapply(current_df,max,na.rm = T)
  #   names(current_max) = paste0(names(current_max),'_max')
  #   
  #   current_mean = sapply(current_df,mean, na.rm = T)
  #   names(current_mean) = paste0(names(current_mean),'_mean')
  #   
  #   current_median = sapply(current_df,median, na.rm = T)
  #   names(current_median) = paste0(names(current_median),'_median')
  #   
  #   current_temp = c(current_max, current_mean, current_median)
  # }
  # 
  # summ_df = append(summ_df,current_temp)









