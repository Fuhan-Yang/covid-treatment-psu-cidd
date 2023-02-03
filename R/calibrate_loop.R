

setwd('..')

source('R/calibrate.R')

smoothed_data = updated_data
for (i in 1:nrow(smoothed_data)) {
	
	smoothed_data$new_cases[i] = mean(updated_data$new_cases[max((i-7), 1):i])
	smoothed_data$new_hosp[i] = mean(updated_data$new_hosp[max((i-7), 1):i])
	smoothed_data$new_deaths[i] = mean(updated_data$new_death[max((i-7), 1):i])
	
}


sae = sum(abs(updated_out$new_hosp[which(updated_out$date >= as.Date('2021-06-01'))] - 
	smoothed_data$new_hosp[which(updated_out$date >= as.Date('2021-06-01'))]), na.rm = TRUE)


num.changes = 0
jump.size = 0
off.tol = 0.1

n.its = 1000

for (it in 1:n.its) {

for (bb in 63:length(beta.strt)) {
		
	beta.new = beta.strt
		
	beta.new[bb] = beta.new[bb] * runif(1, 1-jump.size, 1+jump.size)
	
	test_out = traj.from.params(Z %*% beta.new,
                           params = prms.start,
                           const.params = const.prms.start,
                           non.odesim.params=NULL,
                           introday = introday,
                           tf=end.day,
                           odepath=odepath,
                           loc="RI",
                           symp = NULL,
                           mean.report.rate = mean(rr.daily))

		### get main output ###
		primary_traj = get_ode_main_traj(test_out,scale_up = F)
	
		new_out = data.frame(date = primary_traj$date,
                         new_cases = c(0,diff(primary_traj$cumu_cases,1)),
                         new_hosp = c(0,diff(primary_traj$cumu_hosp,1)),
                         new_deaths = c(0,diff(primary_traj$cumu_death,1)))

	
		sae.new = sum(abs(new_out$new_hosp[which(new_out$date >= as.Date('2021-06-01'))] - 
			smoothed_data$new_hosp[which(updated_out$date >= as.Date('2021-06-01'))]), na.rm = TRUE)
		
		if (sae.new < sae) {
			
			num.changes = num.changes + 1
			sae = sae.new
			beta.strt = beta.new
		
		}
	
}

}

