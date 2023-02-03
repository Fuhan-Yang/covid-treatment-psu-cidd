

#' main function to execute multiple parameter sets
#' 
#' @param start_beta the starting sequence of daily betas
#' @param base_period_start the daynum index of the startday of the baseline period in the start_beta vector
#' @param base_period the number of days considered the baseline period
#' @param extend_period: the number of daily betas to be appended to start_beta
#' @param multiplier: multiplier to the mean betas from the baseline period
#' 
#' @return a sequence of daily betas where the first length(start_beta) elements are start_beta,
#' the next transition_period elements are approximated so that beta values increase/decrease linearly 
#' from the last start_beta to (multiplier * mean(start_beta[-(base_period):])),
#' the last (extend_period - transition_period) elements are multiplier * mean(start_beta[-(base_period):])

# amp is amplitude of seasonality. Represents maximum percentage increase, where 1 is 100% increase
# seas_dur is duration of seasonal increased transmission, in days
# seas_start is the start date, in days of extend_period, of increased transmission

generate_beta2 <- function(start_beta, 
                           base_period_start,
                           base_period = 30, 
                           extend_period = 365, 
                           transition_period = 20,
                           multiplier = 1.0, 
                           amp = 0.10, 
                           seas_dur = 120, 
                           seas_start = 122) {
  
  
  # Because the beta starts at daynum 61, subtract 60 to let the beta period aligns with the actual date #
  base_period_start = base_period_start - 60
  
  if (base_period >= length(start_beta)) {
    stop("base_period must be smaller than the length of start_beta vector", call. = FALSE)
  } else if (transition_period >= extend_period) {
    stop("extend_period must be greater than transition_period", call. = FALSE)
  } else if (seas_start + seas_dur >= 365) {
  	stop("Seasonal increase must start and end within a 365-day period")
  } else if (base_period_start + base_period > length(start_beta)) {
    stop("Sum of base_period_start and base_period must be smaller than the length of start_beta")
  }
 

  ## case start_beta is (n,1) matrix, convert to simple vector
  if ( ! is.null( ncol(start_beta)) ) {
    res <- start_beta[,1]
  } else {
    res <- start_beta
  }

  
  len_start_beta = length(res)
  base_beta = mean(res[base_period_start:(base_period_start + base_period)])
  
  # base_beta = mean(res[(len_start_beta-base_period+1):len_start_beta])
  new_beta = multiplier * base_beta


	# Add annual seasonal multiplier
	seas = rep(1, 365); ind = c(1:365)
	seas[seas_start:(seas_start+seas_dur)] = 1 + 
		amp * sin((2 * pi * (seas_start + seas_dur -ind[seas_start:(seas_start+seas_dur)])) / (2 * seas_dur))
	
	seas = rep(seas, 100); seas = seas[1:(extend_period-1)]
  

  ## transition period
  res <- append(res, approx(c(0,transition_period), c(res[len_start_beta], new_beta), xout=seq(1, transition_period-1))$y  )
  
  

  ## new beta
  res <- append(res, rep(new_beta, extend_period - transition_period))
  
  res[(len_start_beta+1):length(res)] = res[(len_start_beta+1):length(res)] * seas

  ## case start_beta is (n,1) matrix, convert res to (n + extend_period, 1) matrix
  if ( ! is.null( ncol(start_beta)) ){
    dim(res) <- c(length(res), 1)
  }

  return(res)

}

