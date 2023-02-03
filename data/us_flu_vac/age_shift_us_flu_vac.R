rm(list = ls())
library(openxlsx)
library(dplyr)
setwd("~/github/covid19-post-vaccination-burden/")
load('data/us_flu_vac/us_mean_interpolated_flu_vac_across2010_2021_by_age.Rdata')
raw_us_pop = read.xlsx('data/US/us_census/2020_us_census.xlsx')

#### Get the population from the old age groups ###
raw_us_pop = raw_us_pop[-c(1:3),]
raw_us_pop = raw_us_pop %>%
  dplyr::select(`Table.1..Total.U.S..Resident.Population.by.Age,.Sex,.and.Series:.April.1,.2020.(In.thousands)`,
         X4) %>%
  rename(age = `Table.1..Total.U.S..Resident.Population.by.Age,.Sex,.and.Series:.April.1,.2020.(In.thousands)`,
         count = X4) 

us_pop = raw_us_pop

### overall pop is selected regardless of gender ###
us_pop = us_pop[c(1:(which(us_pop$age == 'Male'))-1),]

us_pop$age = as.numeric(us_pop$age) 
us_pop$age = c(0:85)
us_pop$count = as.numeric(us_pop$count) *1000 

### get pop from the vac cov age groups ##
get_pop_by_age = function(age_interval,census = us_pop) {
  
  sum(census$count[which(census$age %in% age_interval)])
  
}
unique(mean_vac_df$age)

### Get the age intervals following the order in the vac cov
vac_age_intervals = list(c(65:69),c(70:85),c(13:17),c(18:19),
                         c(20:49),c(5:12),c(50:64),c(0:4))

pop_by_vac_age = data.frame(age = c('65-69','70-85','13-17','18-19','20-49',
                                    '5-12','50-64','0-4'),
                            pop = NA)

for(i in 1:length(vac_age_intervals)) {
  
  pop_by_vac_age$pop[i] = get_pop_by_age(age_interval = vac_age_intervals[[i]])
                                        
}
pop_by_vac_age

### Get the vac cov after age shifting
all_age_df = data.frame()
for(i in 0:8) {

  this_age_df = data.frame(date = unique(mean_vac_df$weeks),
                           cov = NA,
                           age = paste0('age',i,'0'))
  all_age_df = rbind(all_age_df,this_age_df)

}

vac_by_age = split(all_age_df,all_age_df$age)
names(vac_by_age) = unique(all_age_df$age)

# Get the population by the age groups in the vac data
#### The pop at age20,30,40,50,70,80 would be canceled when calculating the cov so only estimating them
pop_sum_by_old_age = list(age00 = c(pop_by_vac_age$pop[which(pop_by_vac_age$age == '0-4')],
                                    pop_by_vac_age$pop[which(pop_by_vac_age$age == '5-12')] * 5/8),
                          age10 = c(pop_by_vac_age$pop[which(pop_by_vac_age$age == '5-12')] * 3/8,
                                    pop_by_vac_age$pop[which(pop_by_vac_age$age == '13-17')],
                                    pop_by_vac_age$pop[which(pop_by_vac_age$age == '18-19')]),
                          age20 = pop_by_vac_age$pop[which(pop_by_vac_age$age == '20-49')] * 1/3,
                          age30 = pop_by_vac_age$pop[which(pop_by_vac_age$age == '20-49')] * 1/3,
                          age40 = pop_by_vac_age$pop[which(pop_by_vac_age$age == '20-49')] * 1/3,
                          age50 = pop_by_vac_age$pop[which(pop_by_vac_age$age == '50-64')] * 10/15,
                          age60 = c(pop_by_vac_age$pop[which(pop_by_vac_age$age == '50-64')] * 5/15,
                                    pop_by_vac_age$pop[which(pop_by_vac_age$age == '65-69')]),
                          age70 = pop_by_vac_age$pop[which(pop_by_vac_age$age == '70-85')] * 10/16,
                          age80 = pop_by_vac_age$pop[which(pop_by_vac_age$age == '70-85')] * 6/16)


#### get the list by age of the mean vac ###
mean_vac_by_age = split(mean_vac_df, mean_vac_df$age)

### the flu vac cov by old age, assuming the cov is the same across the years in each age group
cov_by_old_age = list(age00 = cbind(mean_vac_by_age$`6 Months - 4 Years`$cov,
                                    mean_vac_by_age$`5-12 Years`$cov),
                      age10 = cbind(mean_vac_by_age$`5-12 Years`$cov,
                                    mean_vac_by_age$`13-17 Years`$cov,
                                    mean_vac_by_age$`18-49 Years`$cov),
                      age20 = mean_vac_by_age$`18-49 Years`$cov,
                      age30 = mean_vac_by_age$`18-49 Years`$cov,
                      age40 = mean_vac_by_age$`18-49 Years`$cov,
                      age50 = mean_vac_by_age$`50-64 Years`$cov,
                      age60 = cbind(mean_vac_by_age$`50-64 Years`$cov,
                                    mean_vac_by_age$`≥65 Years`$cov),
                      age70 = mean_vac_by_age$`≥65 Years`$cov,
                      age80 = mean_vac_by_age$`≥65 Years`$cov)

vac_age_shifter = function(pop_old_age,cov_old_age) {
  
  # vaccinees = pop_age1 * cov_age1 + pop_age2 * cov_age2
  vac = sum(pop_old_age * cov_old_age)
  
  # cov = vaccinees/(pop_age1 + pop_age2)
  vac/sum(pop_old_age)
}

for(age in names(vac_by_age)) {
  
  for(i in 1:nrow(vac_by_age[[age]])) {
    
    if(is.vector(cov_by_old_age[[age]])) {
      vac_by_age[[age]]$cov[i] = vac_age_shifter(pop_old_age = pop_sum_by_old_age[[age]],
                                                 cov_old_age = cov_by_old_age[[age]][i]) 
    } else {
      
      vac_by_age[[age]]$cov[i] = vac_age_shifter(pop_old_age = pop_sum_by_old_age[[age]],
                                                 cov_old_age = cov_by_old_age[[age]][i,]) 
    }
  }
  
}

vac_by_age_df = do.call(rbind,vac_by_age)

# age 20,30,40 have the same lowest flu vac, age 70,80 have the same highest flu vac 
library(ggplot2)
ggplot(vac_by_age_df) + 
  geom_line(aes(x= date, y = cov, color = age))
#ggsave(filename = 'R/data/us_flu_vac/flu_vac_by_age.jpg')

save(vac_by_age_df,file = 'data/us_flu_vac/flu_vac_cov_by_age.Rdata')

# Get the relative coverage through the whole year #
rela_vac_by_age = lapply(1:length(vac_by_age),function(x) {
  
  total_cov = sum(vac_by_age[[x]]$cov)
  rela_cov_df = data.frame(date = vac_by_age[[x]]$date,
                        rela_cov = NA,
                        age = names(vac_by_age)[x])
  for(i in 1:nrow(vac_by_age[[x]])) {
    
    rela_cov_df$rela_cov[i] = vac_by_age[[x]]$cov[i]/total_cov
    
  }
  
  return(rela_cov_df)
})

names(rela_vac_by_age) = names(vac_by_age)
rela_vac_by_age_df = do.call(rbind, rela_vac_by_age)

##### supp fig ######
ggplot(rela_vac_by_age_df) + 
  geom_line(aes(x = date,y = rela_cov,color = age)) + 
  scale_color_discrete(labels = c('0-9','10-19','20-29','30-39','40-49',
                                '50-59','60-69','70-79','>= 80'),
                       name = 'Age Groups') + 
  scale_x_date(date_breaks = '1 month',
               date_labels = '%b',
               expand = c(0,0)) + 
  xlab('') + 
  ylab('Relative Coverage') + 
  theme_bw() + 
  theme(plot.margin = unit(c(0.1,0.1,1.2,0.1),'cm'),
        legend.direction = 'horizontal',
        legend.position = c(0.5,-0.18),
        axis.text = element_text(size = 10)) + 
  guides(color = guide_legend(direction = 'horizontal'))
ggsave(filename = 'data/us_flu_vac/rela_flu_vac_by_age.jpg',
       width = 6,height = 4,units = 'in')

 flu_vac_by_age = rela_vac_by_age

save(flu_vac_by_age,
     file = 'R/data/us_flu_vac/relative_flu_vac_perc_by_age.Rdata')













