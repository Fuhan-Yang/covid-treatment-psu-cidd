# Benefits of near-universal vaccination and treatment access to manage COVID-19 burden in the United States

We evaluated the reduction of COVID-19 burden from vaccination and treatment in a setting with stablized COVID-19 dynamics. We also suggested near-universal vaccination and treatment is the best approach to eliminate COVID-19 through the analysis on the combined strategies.

## Step 1: calibrate the transmission rates to the expected mortality (170,000) in the current season (March 2022 - March 2023)

1. Run `R/calibrate/calibrate_loop.R` to match time-varying transmission parameter to observed hospitalization data
2. Put the calibrated beta into `R/calibrate/baseline_parameters_RI*.R` for three immune durations (365,548, and 730 days). 

## Step 2: generate the annual influenza vaccination trend 

1. Generate the age-specific proportions of weekly counts of influenza vaccines in a year by running `data/us_flu_vac/us_monthly_flu_vac_2009_2021.R` 
2. Shift the age groups to be consistent to our model assuming the probability of getting vaccinated is uniformly distributed across all the ages. This is done by running `data/us_flu_vac/age_shift_us_flu_vac.R`

## Step 3: generate the statewide vaccine and treatment coverage in the current season. 

1. Run `R/functions/calculate_vaccination_coverage.R` 
2. Run `R/functions/get_treat_cov_by_state.R`


## Step 4: generate status quo situation

Run `R/00.status_quo.R` to generate the annual burden under current strategies given specific wane durations, vaccine effectiveness and treatment effectiveness. Note any changes in these parameters require rerun of this script. 

## Step 5: reproduce the results 

### 1. Effectiveness of vaccination: 
 - Run `R/01.vary_vac_cov.R` to generate the total burden between 2025 to 2033 given varying vaccine coverage. 
 - Run `R/summaries_vac_output.R` to generate and plot annual burden between 2025 to 2033 given varying vaccine coverage. 
  
### 2. Effectiveness of treatment: 
 - Run `R/01.vary_treat_cov.R` to generate the total burden between 2025 to 2033 given varying treatment coverage. 
 - Run `R/summaries_treat_output.R` to generate and plot annual burden between 2025 to 2033 given varying treatment coverage. 
  
### 3. Combined strategies: 
 - Run `R/03.vary_vac_treat.R` to generate the total burden between 2025 to 2033 given different combinations of vaccine and treatment coverage. 
 - Run `R/Fig3.R` to generate and plot the annual burden under combined strategies. 
  
### 4. National burden given current statewide strategies
 - Run `04.baseline_by_state.R`
 - Run `R/Fig4.R` to generate and plot the annual burden if current statewide strategies are applied nationally. 
