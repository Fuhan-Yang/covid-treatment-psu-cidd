#### FUNCTIONS TO RUN FUTURE SIMULATIONS ---------------------------------------
library(openxlsx)
#' main function to execute multiple parameter sets
#' 
#' @param change_params data.frame of parameters to change 
#' @param baseline_params list of baseline parameter values including
#' beta, params, const.params, non.odesim.params, introday, tf, and odepath
#' @param output_type character either "summary" (to summarize time series into 
#' public health objectives of interest) or "all" to return the full matrix
#' @param sim_startday, the day where the epidemic burn-in ends
#' @param summary_type choose summary type, details in traj_summary func.
#' @param scale_us logical,if TRUE, scale up to US from RI using monthly ratio of cumu cases
#' between US and RI
#' @param scaled_metrics character lists describing the metrics that will be scaled up to US
#' @return list of length nrow(change_params) 
run_multiple_params <- function(change_params,
                                baseline_params, 
                                output_type = "summary",
                                summary_type = 'max',
                                sim_startday,
                                scale_us = F,
                                scaled_metrics = c('cumu_cases',
                                                   'cumu_hosp',
                                                   'cumu_death')){
  # create output object
  out <- list()
  # use loop for now (can make more efficient if necessary)
  for(i in 1:nrow(change_params)){  
    # update parameters
    params <- param_changer(change_params[i,], baseline_params)
    # run model
    test_out <- traj.from.params(beta = params$beta,
                                 params = params$params,
                                 const.params = params$const.params,
                                 non.odesim.params=params$non.odesim.params,
                                 introday = params$introday,
                                 tf=params$tf,
                                 odepath=params$odepath,
                                 loc=params$loc,
                                 symp = params$symp,
                                 mean.report.rate = params$mean.report.rate)
    # format output
    if(output_type == "primary") {
      if(scale_us) {out[[i]] <- get_ode_main_traj(test_out,
                                                  scale_up = T,
                                                  scaled_metrics = scaled_metrics) 
      } else {out[[i]] <- get_ode_main_traj(test_out)}
    } else if(output_type == "all") {
      out[[i]] <- test_out
    } else if(output_type == 'summary') {
      if(scale_us) {out[[i]] = traj_summary(get_ode_main_traj(test_out,
                                                              scale_up = T,
                                                              scaled_metrics = scaled_metrics),
                                            sim_startday = sim_startday,
                                            options = summary_type)
      } else {out[[i]] = traj_summary(get_ode_main_traj(test_out),
                                      sim_startday = sim_startday,
                                      options = summary_type)}
      
      
    }
  }
  return(out)
}


#' main function to execute multiple parameter sets
#' 
#' @param change_params data.frame of parameters to change 
#' @param baseline_params list of baseline parameter values including
#' beta, params, const.params, non.odesim.params, introday, tf, and odepath
#' @param output_type character either "summary" (to summarize time series into 
#' @param sim_startday, the day where the epidemic burn-in ends
#' public health objectives of interest) or "all" to return the full matrix
#' @param summary_type defines summary type, only available when output_type == 'summary',details in traj_summary
#' @param scale_us logical,if TRUE, scale up to US from RI using monthly ratio of cumu cases
#' between US and RI
#' @param metric the metric will be summarized and scaled 
#' @return list of length nrow(change_params) 
run_multiple_params_age_stratified <- function(change_params,
                                               baseline_params, 
                                               sim_startday,
                                               output_type = "summary",
                                               summary_type = 'max',
                                               scale_us = F,
                                               metric){
  # create output object
  out <- list()
  # use loop for now (can make more efficient if necessary)
  for(i in 1:nrow(change_params)){  
    # update parameters
    params <- param_changer(change_params[i,], baseline_params)
    # run model
    test_out <- traj.from.params(beta = params$beta,
                                 params = params$params,
                                 const.params = params$const.params,
                                 non.odesim.params=params$non.odesim.params,
                                 introday = params$introday,
                                 tf=params$tf,
                                 odepath=params$odepath,
                                 loc=params$loc,
                                 symp = params$symp,
                                 mean.report.rate = params$mean.report.rate)
    # format output
    if(output_type == "primary") {
      if(scale_us) {out[[i]] <- get_age_traj(test_out,
                                             metric = metric,
                                             scale_up = T) 
      } else {out[[i]] <- get_age_traj(test_out,
                                       metric = metric)}
    } else if(output_type == "all") {
      out[[i]] <- test_out
    } else if(output_type == 'summary') {
      if(scale_us) {out[[i]] = traj_summary(get_age_traj(test_out,
                                                         metric = metric,
                                                         scale_up = T),
                                            sim_startday = sim_startday,
                                            options = summary_type)
      } else {out[[i]] = traj_summary(get_age_traj(test_out,
                                                   metric = metric),
                                      sim_startday = sim_startday,
                                      options = summary_type)}
      }
  }
  return(out)
}

#' main function to run single parameter sets combining vaccination and treatment
#' 
#' @param change_params data.frame of parameters(not time-varying) to change
#' @param chagne_vac_cover, annual vaccination coverage, 
#' if NULL, no change in vac cover
#' if a number, fixed annual vaccination coverage for all age groups
#' if a vector, age-varied annual vaccination coverage
#' @inheritParams vac_tv_type
#' @inheritParams vac_cov_ts details in from vac_cover_generator
#' @param baseline_params list of baseline parameter values including
#' beta, params, const.params, non.odesim.params, introday, tf, and odepath
#' @param output_type character either "summary" (to summarize time series into 
#' @param sim_startday, the day where the epidemic burn-in ends
#' public health objectives of interest) or "all" to return the full matrix
#' @param summary_type defines summary type, only available when output_type == 'summary',details in traj_summary
#' @param scale_us logical,if TRUE, scale up to US from RI using monthly ratio of cumu cases
#' between US and RI
#' @param scaled_metrics character lists describing the metrics that will be scaled up to US
#' @return 
#' if output_type == 'primary', return a data.frame of trajs of cumu cases,hosp, and deaths, can be scaled to US
#' if output_type == 'summary', return a vector of summarized trajs of cumu cases, hosp, and deaths, can be scaled to US
#' if output_type == 'all', return the raw output from cpp, cannot be scaled 
run_single_param_comb = function(change_params,
                                 change_vac_cover = NULL,
                                 pop_by_age,
                                 vac_tv_type = 'constant_increase',
                                 vac_cov_ts = NULL,
                                 baseline_params,
                                 output_type = 'summary',
                                 summary_type = 'max',
                                 sim_startday,
                                 scale_us = F,
                                 scaled_metrics = c('cumu_cases',
                                                    'cumu_hosp',
                                                    'cumu_death')) {
  out = list()
  ### update one-element parameters
  params = param_changer(change_params,baseline_params)
  ### update time-varying vaccine coverage
  ## constant vac cover across ages
  if(is.null(change_vac_cover)) {
    params = params
  } else if(length(change_vac_cover == 1)) {
    params = vac_cover_generator(cover = change_vac_cover,
                                 pop_by_age = pop_by_age,
                                 age_cov_varied = F,
                                 vac_tv_type = vac_tv_type,
                                 vac_cov_ts = vac_cov_ts,
                                 baseline_params = params)
    
  } else if(length(change_vac_cover > 1)) {
    params = vac_cover_generator(cover = change_vac_cover,
                                 pop_by_age = pop_by_age,
                                 age_cov_varied = T,
                                 cover_by_age = change_vac_cover, 
                                 vac_tv_type = vac_tv_type,
                                 vac_cov_ts = vac_cov_ts,
                                 baseline_params = params)
    
  }
  # run model
  test_out <- traj.from.params(beta = params$beta,
                               params = params$params,
                               const.params = params$const.params,
                               non.odesim.params=params$non.odesim.params,
                               introday = params$introday,
                               tf=params$tf,
                               odepath=params$odepath,
                               loc=params$loc,
                               symp = params$symp,
                               mean.report.rate = params$mean.report.rate)
  # format output
  if(output_type == "primary"){
    if(scale_us) {
      out = get_ode_main_traj(test_out,
                              scale_up = T,
                              scaled_metrics = scaled_metrics)
    } else {
      out <- get_ode_main_traj(test_out)  
    }
  } else if(output_type == "all"){
    out <- test_out
  } else if(output_type == 'summary') {
    if(scale_us) {
      out = traj_summary(get_ode_main_traj(test_out,
                                           scale_up = T,
                                           scaled_metrics = scaled_metrics),
                         sim_startday = sim_startday,
                         options = summary_type)
    } else {
      
      out = traj_summary(get_ode_main_traj(test_out),
                         sim_startday = sim_startday,
                         options = summary_type)
    }
  }
  
  
  return(out)
}


#' main function to run single parameter sets combining vaccination and treatment, return age-specific traj
#' 
#' @param change_params data.frame of parameters(not time-varying) to change
#' @param chagne_vac_cover, annual vaccination coverage, 
#' if NULL, no change in vac cover
#' if a number, fixed annual vaccination coverage for all age groups
#' if a vector, age-varied annual vaccination coverage
#' @inheritParams vac_startday
#' @inheritParams vac_endday
#' @inheritParams vac_tv_type
#' @inheritParams vac_cov_ts details in from vac_cover_generator
#' @param baseline_params list of baseline parameter values including
#' beta, params, const.params, non.odesim.params, introday, tf, and odepath
#' @param output_type character either "summary" (to summarize time series into 
#' public health objectives of interest) or "all" to return the full matrix
#' @param sim_startday, the day where the epidemic burn-in ends
#' @param summary_type defines summary type, only available when output_type == 'summary',details in traj_summary
#' @param scale_us logical,if TRUE, scale up to US from RI using monthly ratio of cumu cases
#' between US and RI
#' @param metrics traj metrics; choose from 'cases','hosp','deaths'
#' @return 
#' if output_type == 'primary', return a data.frame of trajs of cumu cases,hosp, and deaths, can be scaled to US
#' if output_type == 'summary', return a vector of summarized trajs of cumu cases, hosp, and deaths, can be scaled to US
#' if output_type == 'all', return the raw output from cpp, cannot be scaled 
run_single_param_comb_age_stratified = function(change_params,
                                 change_vac_cover = NULL,
                                 pop_by_age,
                                 vac_startday, 
                                 vac_endday,
                                 sim_startday,
                                 vac_tv_type = 'constant_increase',
                                 vac_cov_ts = NULL,
                                 baseline_params,
                                 output_type = 'summary',
                                 summary_type = 'max',
                                 scale_us = F,
                                 metric) {
  out = list()
  ### update one-element parameters
  params = param_changer(change_params,baseline_params)
  ### update time-varying vaccine coverage
  ## constant vac cover across ages
  if(is.null(change_vac_cover)) {
    params = params
  } else if(length(change_vac_cover == 1)) {
    params = vac_cover_generator(cover = change_vac_cover,
                                 pop_by_age = pop_by_age,
                                 age_cov_varied = F,
                                 vac_tv_type = vac_tv_type,
                                 vac_cov_ts = vac_cov_ts,
                                 baseline_params = params)
    
  } else if(length(change_vac_cover > 1)) {
    params = vac_cover_generator(cover = change_vac_cover,
                                 pop_by_age = pop_by_age,
                                 age_cov_varied = T,
                                 cover_by_age = change_vac_cover, 
                                 vac_tv_type = vac_tv_type,
                                 vac_cov_ts = vac_cov_ts,
                                 baseline_params = params)
    
  }
  # run model
  test_out <- traj.from.params(beta = params$beta,
                               params = params$params,
                               const.params = params$const.params,
                               non.odesim.params=params$non.odesim.params,
                               introday = params$introday,
                               tf=params$tf,
                               odepath=params$odepath,
                               loc=params$loc,
                               symp = params$symp,
                               mean.report.rate = params$mean.report.rate)
  # format output
  if(output_type == "primary"){
    if(scale_us) {
      out = get_age_traj(test_out,
                         metric = metric,
                         scale_up = T)
    } else {
      out <- get_age_traj(test_out,
                          metric = metric,
                          scale_up = F)
    }
  } else if(output_type == "all"){
    out <- test_out
  } else if(output_type == 'summary') {
    if(scale_us) {
      out = traj_summary(get_age_traj(test_out,
                                      metric = metric,
                                      scale_up = T),
                         sim_startday = sim_startday,
                         options = summary_type)
    } else {
      
      out = traj_summary(get_age_traj(test_out,
                                      metric = metric,
                                      scale_up = F),
                         sim_startday = sim_startday,
                         options = summary_type)
    }
  }
  
  return(out)
}





#' helper function 'param_changer'
#' for now, only change the parameters in params and const.params
#' 
#' @param change_params a named vector for multiple parameters 
#' @param baseline_params previously estimated params plus future beta, 
#' differ by state
#' @return parameter list with change_param values updated from baseline_params
param_changer <- function(change_params, baseline_params) {
  new_params <- baseline_params
  # check the length of change_params
  
  # check for parameters to change in baseline_params$params
  p_loc <- match(names(change_params), names(new_params$params))
  new_params[["params"]][p_loc[!is.na(p_loc)]] <- unlist(change_params[!is.na(p_loc)])
  # and repeat for baseline_params$const.params
  p_loc <- match(names(change_params), names(new_params$const.params))
  new_params[["const.params"]][p_loc[!is.na(p_loc)]] <- unlist(change_params[!is.na(p_loc)]) 
  
  return(new_params)
}


#' helper function of 'vac_param_append': convert vec to string
#' @vec should be the input of new vac coverage/vac endday
#' @return should be string that can be read by C++
vec_to_string = function(vec) {
  
  vec = data.frame(vec)
  # no delimiter is '\t\, coerce to '\n\' #
  string = readr::format_delim(vec,delim = '\n',col_names = F)
  # replace '\n' to '\t' #
  string = gsub('\n','\t',string)
  # remove the last extra '\t' #
  string = gsub('(.*)\t','\\1',string)
  
  return(string)
}


#### Update needed: differ from state and national pop ####
#' @param loc specifies the location of population, only US and RI is available
#' read 2019 census file downloaded from https://www.census.gov/data/tables/time-series/demo/popest/2010s-state-detail.html
#' shift the age to the desired age groups
#' @return data frame of population by age 
get_pop_by_age = function(loc) {
  
  prev_dir = getwd()
  if(prev_dir != "/Users/fuhanyang/github/covid19-post-vaccination-burden") {
    
    setwd("/Users/fuhanyang/github/covid19-post-vaccination-burden")
  }
  if(loc == 'RI') {
    
    census = read.xlsx('data/RI/RI-census-by-age.xlsx',
                       startRow = 3)
    
    census2019 = data.frame(cbind(census[,1],census[,35]))
    census2019 = census2019[3:nrow(census2019),]
    census2019[,2] = as.numeric(census2019[,2])
    
    all_needed_age_names = c('.Under 5 years','.5 to 9 years','.10 to 14 years',
                             '.15 to 19 years','.20 to 24 years','.25 to 29 years',
                             '.30 to 34 years','.35 to 39 years','.40 to 44 years',
                             '.45 to 49 years','.50 to 54 years','.55 to 59 years',
                             '.60 to 64 years','.65 to 69 years','.70 to 74 years',
                             '.75 to 79 years','.80 to 84 years','.85 years and over')
    pop2019 = unique(census2019[which(census2019[,1] %in% all_needed_age_names),])
    
    pop_by_age = data.frame(
      age0 = NA,
      age10 = NA,
      age20 = NA,
      age30 = NA,
      age40 = NA,
      age50 = NA, 
      age60 = NA,
      age70 = NA,
      age80 = NA)                   
    
    pop_by_age$age0 = sum(pop2019[pop2019[,1] == '.Under 5 years',2],
                          pop2019[pop2019[,1] == '.5 to 9 years',2])
    pop_by_age$age10 = sum(pop2019[pop2019[,1] == '.10 to 14 years',2],
                           pop2019[pop2019[,1] == '.15 to 19 years',2])
    pop_by_age$age20 = sum(pop2019[pop2019[,1] == '.20 to 24 years',2],
                           pop2019[pop2019[,1] == '.25 to 29 years',2])
    pop_by_age$age30 = sum(pop2019[pop2019[,1] == '.30 to 34 years',2],
                           pop2019[pop2019[,1] == '.35 to 39 years',2])
    pop_by_age$age40 = sum(pop2019[pop2019[,1] == '.40 to 44 years',2],
                           pop2019[pop2019[,1] == '.45 to 49 years',2])
    pop_by_age$age50 = sum(pop2019[pop2019[,1] == '.50 to 54 years',2],
                           pop2019[pop2019[,1] == '.55 to 59 years',2])
    pop_by_age$age60 = sum(pop2019[pop2019[,1] == '.60 to 64 years',2],
                           pop2019[pop2019[,1] == '.65 to 69 years',2])
    pop_by_age$age70 = sum(pop2019[pop2019[,1] == '.70 to 74 years',2],
                           pop2019[pop2019[,1] == '.75 to 79 years',2])
    pop_by_age$age80 = sum(pop2019[pop2019[,1] == '.80 to 84 years',2],
                           pop2019[pop2019[,1] == '.85 years and over',2])
  } else if (loc == 'US') {
    
    census = read.csv('data/US/us_census/US_population_2020_estimates.csv')
    pop_by_age = as.numeric(t(census[,2]))
    names(pop_by_age) = census[,1]
    
  } else {
    
    stop('loc not available')
  }
  setwd(prev_dir)
  return(pop_by_age)
}


#' internal function of 'vac_cover'
#' @change_params the updated vac params, should be a df with vac_endday and at least 1 vac_coverage
#' @baseline_params the previously estimated params, differed by state
#' @return the updated params comb, ready for traj.from.params
vac_param_append = function(change_params, baseline_params) {
  
  new_params <- baseline_params
  # when change_params are longer than 1, it should be in vaccine coverage data
  # vac cover params should include vac-endday and vaccinees, the ncol should be >1
  vac_cover_name = c('tv-vaccinees-00','tv-vaccinees-10','tv-vaccinees-20',
                     'tv-vaccinees-30','tv-vaccinees-40','tv-vaccinees-50',
                     'tv-vaccinees-60','tv-vaccinees-70','tv-vaccinees-80')
  # vac change params should have at least 2 same colnames as new_params
  if(!any(colnames(change_params) == 'tv-vaccinees-endday') |
     !any(colnames(change_params) %in% vac_cover_name)) {
    stop('change_params name is inconsistent with vac.params names, check again.')
  }
  
  # only change the vac cover after the time when vac data is available
  baseline_vac_enddays = as.numeric(unlist(strsplit(new_params$const.params['tv-vaccinees-endday'],
                                                    split = '\t')))
  if(min(change_params['tv-vaccinees-endday']) <= max(baseline_vac_enddays)) {
    stop('the endday of vac in change_params should be later than vac data')
  }
  
  for(i in 1:ncol(change_params)) { 
    
    loc = which(names(new_params$const.params) == colnames(change_params)[i])
    new_params$const.params[loc] = paste0(new_params[['const.params']][loc],'\t',
                                          vec_to_string(change_params[,i]))
  }
  
  return(new_params)
  
}

#### pop by age from census is needed ####
#' create tv-vac-cover based on total coverage
#' @cover the percentage of vaccinated population in the total population
#' @startday the day when vaccine campaign starts: this should be the firstday of the week after the available vaccine data
#' @endday the day when vaccine campaign ends
#' @age_cov_varied if the vac coverage differs by age? 
#' @vac_tv_type the time-varying pattern of vac coverage, assuming it is constantly increasing or other scenarios
#' @vac_cov_ts the vaccine coverage time series simulating other vac coverage when vac_tv_type == 'scenarios',should be a list with 9 elements
#' @baseline_params the previous estimates, required in traj.from.params
#' @return_vac return the simulated vac cov ts, only for debug
#' @return the params for traj.params with updated vac coverage
vac_cover_generator = function(
    cover,
    pop_by_age,
    age_cov_varied = F,
    cover_by_age = NULL,
    vac_tv_type = c('constant_increase','scenarios'),
    vac_cov_ts = NULL, 
    return_vac = F,
    baseline_params) {
  
  ## if number of total vaccinees differs by age ##
  if(age_cov_varied == T) {
    if(is.null(cover_by_age)) {stop('age-specific coverage is needed if coverage is varied by age')}
    ### need to change this when age 0-9 is included 
    if(length(cover_by_age) !=9) {stop('There should be 9 coverages for 0-80 age groups.')}
    if(length(unique(cover_by_age)) ==1) {stop('coverage is the same across the ages, specify cover instead')}
    
    ### cover_by_age must have the same age order as pop_by_age ###
    total_vac_by_age = round(pop_by_age*cover_by_age)
    
  } else {
    ### same coverage across all the age groups
    total_vac_by_age = round(pop_by_age*cover)
  }
  total_vac_by_age = as.numeric(total_vac_by_age)
  
  ### need to change this when age 0-9 is included 
  names(total_vac_by_age) = paste0('age',seq(0,80,length.out = 9))
  
  ### vac cov starts to generate the next week of available vac data, and stops at the end of simulation ###
  startday = max(as.numeric(unlist(strsplit(baseline_params$const.params[["tv-vaccinees-endday"]],'\t')))) + 7
  endday = baseline_params$tf
  
  # this vaccine pattern assumes the population get vaccinated constantly through the whole time series with no age differences #
  if(vac_tv_type == 'constant_increase') {
    
    ## time-varying vac cover ##
    ## vac data is reported weekly ##
    tv_vac_by_age = data.frame(time = seq(startday, endday, 7))
    
    tv_vac_by_age$wts = 1/nrow(tv_vac_by_age)
    vac_by_age_ts = lapply(1:length(total_vac_by_age),function(x) {round(tv_vac_by_age$wts*total_vac_by_age[x])})
    vac_by_age_ts = do.call(cbind, vac_by_age_ts)
    tv_vac_by_age$wts = NULL
  }
  
  # this vaccine pattern assumes different realistic vaccine scenarios, informed by vac_cov_ts #
  if(vac_tv_type == 'scenarios') {
    
    ### age groups should be ordered increasingly ###
    if(length(vac_cov_ts) != 9) {stop('There should be 9 timeseries for 9 age groups.')} 
    
    vac_by_age_ts = lapply(1:length(vac_cov_ts),function(x) {
      
      vac_cov_ts[[x]]$age = NULL
      vac_cov_ts[[x]]$week_ind = findInterval(vac_cov_ts[[x]]$date,vac_cov_ts[[x]]$date)
      
      # convert daynum to date #
      startdate = as.Date('2020-01-01') + startday
      enddate = as.Date('2020-01-01') + endday
      
      ## find which week index does the startdate belong to ## 
      fake_startdate = as.Date(paste0('2020-',format(startdate,'%m-%d')))
      startweek_ind = findInterval(fake_startdate,vac_cov_ts[[x]]$date)
      
      ## find which week index does the enddate belong to ## 
      fake_enddate = as.Date(paste0('2020-',format(enddate,'%m-%d')))
      endweek_ind = findInterval(fake_enddate,vac_cov_ts[[x]]$date)
      
      
      ### the sim is start_remains in the first sim year + years_in_between + end_remains in the last year ###
      
      ##  The vac cov is the same for all days within one week, so we use weekly index to define vac cov ##
      start_remains = vac_cov_ts[[x]][startweek_ind:nrow(vac_cov_ts[[x]]),]
      start_remains$date = seq(startdate,as.Date(startdate + 7 * (nrow(vac_cov_ts[[x]]) - startweek_ind)),'7 days')
      
      end_remains = vac_cov_ts[[x]][1:endweek_ind,]
      end_remains$date = seq(as.Date(enddate - 7 * (endweek_ind - 1)),enddate,'7 days')
      
      ### When there is at least 1 year between the startdate and enddate ###
      if(max(start_remains$date) + 7 < min(end_remains$date) - 7) {
      years_inbetween = data.frame(date = seq(max(start_remains$date) + 7, min(end_remains$date) - 7,'7 days'),
                                   week_ind = NA, 
                                   rela_cov = NA)
      
      years_inbetween$fake_date = as.Date(paste0('2020-',format(years_inbetween$date,'%m-%d')))
      years_inbetween$week_ind = findInterval(years_inbetween$fake_date,vac_cov_ts[[x]]$date)
      
      for(i in 1:nrow(years_inbetween)) {
        
        years_inbetween$rela_cov[i] = vac_cov_ts[[x]]$rela_cov[which(vac_cov_ts[[x]]$week_ind == years_inbetween$week_ind[i])]
      }
      
      years_inbetween$fake_date = NULL
      
      whole_sim = rbind(start_remains,years_inbetween,end_remains)
      
      ### When the interval between the startdate and enddate is no more than 1 year ###
      } else {
        
        whole_sim = rbind(start_remains,end_remains)
      }
      
      whole_sim$daynum = as.numeric(difftime(whole_sim$date,as.Date('2020-01-01')))
      whole_sim$week_ind = NULL
      whole_sim$date = NULL
      
      return(whole_sim)
      
    })
    
    tv_vac_by_age = data.frame(time = vac_by_age_ts[[1]]$daynum)
    vac_by_age_ts = lapply(1:length(vac_by_age_ts),function(x){round(vac_by_age_ts[[x]]$rela_cov*total_vac_by_age[x])})
    vac_by_age_ts = do.call(cbind,vac_by_age_ts)
    
  }
  colnames(vac_by_age_ts) = names(total_vac_by_age)
  
  tv_vac_by_age = data.frame(tv_vac_by_age,vac_by_age_ts)
  
  colnames(tv_vac_by_age) = c('tv-vaccinees-endday','tv-vaccinees-00',
                              'tv-vaccinees-10','tv-vaccinees-20',
                              'tv-vaccinees-30','tv-vaccinees-40','tv-vaccinees-50',
                              'tv-vaccinees-60','tv-vaccinees-70','tv-vaccinees-80')
  
  if(return_vac == T) 
    
  {return(tv_vac_by_age)} 
  
  else {
    
    updated_vac_params = vac_param_append(change_params = tv_vac_by_age,
                                          baseline_params = baseline_params)
    return(updated_vac_params)
  }
  
  
}





