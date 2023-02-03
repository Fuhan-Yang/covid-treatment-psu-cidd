library(openxlsx)
library(dplyr)
setwd("~/github/covid19-post-vaccination-burden")
raw_us_pop = read.xlsx('data/US/us_census/2020_us_census.xlsx')

#### Get the population from the old age groups ###
raw_us_pop = raw_us_pop[-c(1:3),]
raw_us_pop = raw_us_pop %>%
  select(`Table.1..Total.U.S..Resident.Population.by.Age,.Sex,.and.Series:.April.1,.2020.(In.thousands)`,
         X4) %>%
  rename(age = `Table.1..Total.U.S..Resident.Population.by.Age,.Sex,.and.Series:.April.1,.2020.(In.thousands)`,
         count = X4) 

us_pop = raw_us_pop

### overall pop is selected regardless of gender ###
us_pop = us_pop[c(1:(which(us_pop$age == 'Male'))-1),]

us_pop$age = as.numeric(us_pop$age) 
us_pop$age = c(0:85)
us_pop$count = as.numeric(us_pop$count) *1000 

age_shifted_groups = list(0:9,10:19,20:29,30:39,40:49,50:59,60:69,70:79,80:85)

age_shifted_pop = c()

for(i in 1:length(age_shifted_groups)) {
  
  age_shifted_pop[i] = sum(us_pop$count[which(us_pop$age %in% age_shifted_groups[[i]])])
}
names(age_shifted_pop) = paste0('age',seq(0,80,10))

# age_shifted_pop = data.frame(age_groups = paste0('age',seq(0,80,10)),
#                              pop = NA)
# for(i in 1:nrow(age_shifted_pop)) {
#   
#   age_shifted_pop$pop[i] = sum(us_pop$count[which(us_pop$age %in% age_shifted_groups[[i]])])
#   
# }

write.csv(age_shifted_pop,file = 'data/US/us_census/US_population_2020_estimates.csv')
