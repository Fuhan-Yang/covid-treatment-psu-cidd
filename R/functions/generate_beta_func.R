# beta <- Z %*% beta.strt

#' main function to execute multiple parameter sets
#' 
#' @param start_beta the starting sequence of daily betas
#' @param base_period the number of days considered the baseline period, counting backwards from the last day of the start_beta
#' @param extend_period: the number of daily betas to be appended to start_beta
#' @param multiplier: multiplier to the mean betas from the baseline period
#' 
#' @return a sequence of daily betas where the first length(start_beta) elements are start_beta,
#' the next transition_period elements are approximated so that beta values increase/decrease linearly 
#' from the last start_beta to (multiplier * mean(start_beta[-(base_period):])),
#' the last (extend_period - transition_period) elements are multiplier * mean(start_beta[-(base_period):])
generate_beta <- function(start_beta, 
                          base_period = 30, 
                          extend_period = 365, 
                          transition_period = 20,
                          multiplier = 1.0){
  
  if (base_period >= length(start_beta)) {
    stop("base_period must be smaller than than the length of start_beta vector", call. = FALSE)
  } else if (transition_period >= extend_period) {
    stop("extend_period must be greater than transition_period", call. = FALSE)
  }
  
  ## case start_beta is (n,1) matrix, convert to simple vector
  if ( ! is.null( ncol(start_beta)) ){
    res <- start_beta[,1]
  } else {
    res <- start_beta
  }
  
  len_start_beta = length(res)
  base_beta = mean(res[(len_start_beta-base_period+1):len_start_beta])
  new_beta = multiplier * base_beta
  
  ## transition period
  res <- append(res, approx(c(0,transition_period), c(res[len_start_beta], new_beta), xout=seq(1, transition_period-1))$y  )
  
  ## new beta
  res <- append(res, rep(new_beta, extend_period - transition_period))
  
  ## case start_beta is (n,1) matrix, convert res to (n + extend_period, 1) matrix
  if ( ! is.null( ncol(start_beta)) ){
    dim(res) <- c(length(res), 1)
  }
  
  return(res)
}
