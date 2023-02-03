library(dplyr)
library(data.table)
library(lubridate)
library(ggplot2)
library(scales)



#### LOAD OBSERVED VACCINATION DATA --------------------------------------------
# data pulled from CDC COVID-19 vaccination database
# https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-Jurisdi/unsk-b7fc
# ^ visit this URL for metadata ^
cdc_path <- "https://data.cdc.gov/api/views/unsk-b7fc/rows.csv?accessType=DOWNLOAD" 

# CDC vaccination database
vacc_raw <- fread(cdc_path)
vacc_raw[, Date := mdy(Date)]

#### CALCULATE VACCINATION COVERAGE --------------------------------------------

# set cutoff dates
# use Dec 1, 2021 - November 30, 2022
first_date = as.Date("2021-12-01")
last_date = as.Date("2022-11-30")

## calculating % of individuals in US that have 
## received one dose in the past year (between first_date and last_date)
## CDC fields do not clearly subset so pull a few key metrics during this
## time period, and use these to estimate conservative and optimistic coverage 
## Note: not using age specific fields, because they are not provided consistently
## for all of the fields we need

# 1. total administered doses (i.e., shots in arms)
admin_full = vacc_raw %>%
  filter(Location == "US", Date == last_date) %>%
  pull(Administered) - 
  vacc_raw %>% 
  filter(Location == "US", Date == first_date) %>%
  pull(Administered)
admin_full/1E6
# 2. total first doses administered (i.e., first dose in primary series) 
first_doses = vacc_raw %>%
  filter(Location == "US", Date == last_date) %>%
  pull(Administered_Dose1_Recip) - 
  vacc_raw %>% 
  filter(Location == "US", Date == first_date) %>%
  pull(Administered_Dose1_Recip)
first_doses/1E6
# 3. number of individuals that have completed two-course series 
series_complete = vacc_raw %>%
  filter(Location == "US", Date == last_date) %>%
  pull(Series_Complete_Yes) - 
  vacc_raw %>% 
  filter(Location == "US", Date == first_date) %>%
  pull(Series_Complete_Yes)
series_complete/1E6
# 4. total additional doses administered
additional_admin_complete = vacc_raw %>%
  filter(Location == "US", Date == last_date) %>%
  pull(Additional_Doses) - 
  vacc_raw %>% 
  filter(Location == "US", Date == first_date) %>%
  pull(Additional_Doses)
additional_admin_complete/1E6
# 5. second boosters (non-bivalent)
second_boost = vacc_raw %>%
  filter(Location == "US", Date == last_date) %>%
  pull(Second_Booster) - 
  vacc_raw %>% 
  filter(Location == "US", Date == first_date) %>%
  # second booster data starts after first_date
  mutate(Second_Booster = ifelse(is.na(Second_Booster), 
                                 0, Second_Booster)) %>%
  pull(Second_Booster)
second_boost/1E6
# 6. bivalent boosters
bival_boost = vacc_raw %>%
  filter(Location == "US", Date == last_date) %>%
  pull(Administered_Bivalent) - 
  vacc_raw %>% 
  filter(Location == "US", Date == first_date) %>%
  # second booster data starts after first_date
  mutate(Administered_Bivalent = ifelse(is.na(Administered_Bivalent), 
                                        0, Administered_Bivalent)) %>%
  pull(Administered_Bivalent)
bival_boost/1E6

## using these values, calculate coverage estimates
# US population size
# use two fields in vacc_raw to get this
US_pop = vacc_raw[Date == last_date & Location == "US", 
                  .(Series_Complete_Yes, Series_Complete_Pop_Pct)] %>%
  .[, pop_size := Series_Complete_Yes/(Series_Complete_Pop_Pct/100)] %>% pull(pop_size)

# optimistic coverage estimate
# assumes the only duplicate doses were individuals completing two-course series
opt_cov = (admin_full - series_complete)/US_pop
round(opt_cov,3)*100

# conservative coverage estimate
# assumes that individuals completing two-course series and individuals recieveing
# a second booster all had received another dose at some point during the time
# period (i.e., these individuals are double counted)
cons_cov = (admin_full - series_complete - second_boost)/US_pop
round(cons_cov,3)*100


# a plot to demonstrate the two values we're working with here (for the US)
vacc_raw[Location == "US", 
         .(Date, Administered, Series_Complete_Yes, Administered_Dose1_Recip, 
           Additional_Doses, Second_Booster, Administered_Bivalent)] %>%
  melt(c("Date")) %>%
  mutate(variable = factor(variable, 
                           levels = c("Administered", "Administered_Dose1_Recip", 
                                      "Series_Complete_Yes", "Additional_Doses", 
                                      "Second_Booster", "Administered_Bivalent"))) %>%
ggplot(aes(x = Date)) + 
  geom_line(aes(y = value, color = variable), size = 1) + 
  geom_vline(xintercept = first_date, linetype = "dashed") + 
  geom_vline(xintercept = last_date, linetype = "dashed") + 
  guides(color=guide_legend(nrow=3)) +
  scale_color_discrete(labels = c("Administered", 
                                  "First dose", 
                                  "Completed primary series", 
                                  "Additional booster dose (non-bivalent)", 
                                  "Second booster (non-bivalent)", 
                                  "Bivalent booster")) +
  scale_x_date(date_labels = "%B %Y") +
  scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6), 
                     name = "US vaccine doses") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank())
# ggsave("figs/suppFig.vaccination_uptake.pdf", width = 6, height = 6.5)

#### IMPLEMENT FOR ALL LOCATIONS -----------------------------------------------
# need to calculate state-level population sizes
# hack this by dividing two fields in the vacc data
state_pop <- vacc_raw[Date == last_date, .(Location, Series_Complete_Yes, Series_Complete_Pop_Pct)] %>%
  .[, pop_size := Series_Complete_Yes/(Series_Complete_Pop_Pct/100) , by = .(Location) ]


pct_vacc_by_state <- vacc_raw[Date %in% as.Date(c(first_date, last_date)), 
                             .(Date, Location, Administered, Series_Complete_Yes, 
                               Administered_Dose1_Recip, Additional_Doses, 
                               Second_Booster, Administered_Bivalent)] %>%
  melt(c("Date", "Location")) %>%
  # some fields do not have data on start_date, so use 0 
  .[, value := ifelse(is.na(value), 0, value)] %>%
  # subset start_date and end_date into two discrete fields 
  .[, which_date := ifelse(Date == first_date, "first_date", "end_date")] %>%
  data.table::dcast(variable + Location ~ which_date, value.var = "value") %>%
  .[, difference := end_date - first_date] %>% 
  # recast
  data.table::dcast(Location ~ variable, value.var = "difference") %>%
  # optimistic = admin_full - series_complete
  # note: second booster (and other info) does not seem to be available at state level
  # so we only calculate optimistic here
  .[, opt := Administered - Series_Complete_Yes] %>%
  # add in population sizes
  .[state_pop[, .(Location, pop_size)], on = .(Location)] %>%
  .[, ":=" (opt_pct_vacc = opt/pop_size)] %>%
  # remove non-state Locations
  .[!(Location %in% c("BP2", "DD2", "VA2", "IH2", "PR", "FM", "MH", "AS", "GU", "PW", "VI"))]
save(pct_vacc_by_state,file = 'R/cpp-v10/Rdata/vac_cov_baseline_by_state.Rdata')


# plot to see the coverage values by state
pct_vacc_by_state[, .(Location, opt_pct_vacc, pop_size)] %>%
  melt(c("Location", "pop_size")) %>%
  ggplot(aes(x = reorder(Location, -value), y = value * 100)) + 
  geom_point(size = 3) +
  # geom_text(aes(label = paste0(round(value,2)*100, "%")),
  #           size = 2.5, color = "white") +
  scale_y_continuous(n.breaks = 10) + 
  xlab('State') + 
  ylab('Vaccination Coverage (%)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = 'figs/supp_fig_vac_cov_by_state.jpg',width = 8,height = 4,
       units = 'in')



#### TRY DOING THE US LEVEL ESTIMATE BY VACC % ####
# using: covid19-post-vaccination-burden/R/Rdata/vac_treat_summaries.Rdata

# estimate cumualtive deaths for each state based on vaccination coverage
# using linear interpolation between the values we obtain from simulation
treat_covs_to_plot <- c(0.2, 0.5, 0.8)
deaths_by_state <- list()
for(i in 1:length(treat_covs_to_plot)){
  a <- all_param_combs_list[["neutral"]] %>%
    filter(treat == treat_covs_to_plot[i])
  a <- approx(a$vac, a$cumu_death, pct_vacc_by_state$pct_vacc)
  deaths_by_state[[i]] <- data.frame(Location = pct_vacc_by_state$Location, 
                                     pct_vacc = pct_vacc_by_state$pct_vacc, 
                                     est_death = a$y, 
                                     treat_cov = treat_covs_to_plot[i])
}
deaths_by_state <- data.table::rbindlist(deaths_by_state)

## US roll up estimates
deaths_by_state %>%
  group_by(treat_cov) %>%
  summarize(us_deaths = sum(est_death),
            annual_US_deaths_M = sum(est_death)/1000)

# plot to illustrate
states_to_plot <- c("VI", "TX", "PA", "NY")
vacc_points <- deaths_by_state %>% 
  filter(Location %in% states_to_plot)

x_labs <- data.frame(val = c(seq(0,1, 0.25), vacc_points$pct_vacc %>% unique()), 
                     lab = c(paste0(seq(0,1, 0.25)*100,"%"), as.character(vacc_points$Location %>% unique())))
x_labs <- x_labs[order(x_labs$val), ]

all_param_combs_list[["neutral"]] %>%
  filter(treat %in% c(0.2, 0.5, 0.8)) %>%
  ggplot(aes(x = vac, y = cumu_death/10, color = as.factor(treat))) + 
  geom_line(size = 1) + 
  geom_point(data = vacc_points, 
             aes(x = pct_vacc, y = est_death/10, color = as.factor(treat_cov)), 
             size = 3) +
  geom_text(data = vacc_points, 
             aes(x = pct_vacc, y = est_death/10, label = Location), 
             size = 1.5, color = "white") +
  scale_color_manual(values =RColorBrewer::brewer.pal(4, "Purples")[2:4], 
                     labels = c("20%", "50%", "80%"),
                     name = "treatment coverage") +
  scale_x_continuous(expand = c(0,0), 
                     labels = percent,
                     name = "vaccination coverage") +
  scale_y_continuous(expand = c(0,0),
                     labels = label_number(suffix = "K", scale = 1e-3), 
                     name = "average yearly deaths") +
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.grid.minor = element_blank())
ggsave("sample_state_rollup.pdf", width = 4, height = 4)


