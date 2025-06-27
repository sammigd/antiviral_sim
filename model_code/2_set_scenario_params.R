

##################################################
# CALC BETAS #
##################################################

#SAR HARDCODED IN FROM UNPUBLISHED SAR IMPUTATION WORK
#for age groups 0-4, 5-17, 18-49, 50-64,50+
SAR_flutes = c(runif(1, 0.11, 0.34), # ranges from sinead appendix table
               runif(1, 0.12, 0.44),
               runif(1, 0.06, 0.14), 
               runif(1, 0.07, 0.18),
               runif(1, 0.02, 0.09))

SAR_flutes = c(runif(1, 0.21 - 0.01, 0.21 + 0.01), # flutes est + arbitrary ranges
               runif(1, 0.14 - 0.01, 0.14 + 0.01),
               runif(1, 0.19 - 0.01, 0.19 + 0.01), 
               runif(1, 0.25 - 0.01, 0.25 + 0.01),
               runif(1, 0.25 - 0.01, 0.25 + 0.01))

SAR_flutes = c(0.21, 0.14, 0.19, 0.25, 0.25)

#########################################
scenario_tab = read.csv('inputs/scenarios.csv', row.names = 1)
epi_scenarios = names(scenario_tab)

av_transmissibility_adj   = scenario_tab['av_transmissibility_adj', epi_scenario]
sar_adj_noav              = scenario_tab['sar_adj_noav', epi_scenario]
p_adherence               = scenario_tab['p_adherence', epi_scenario]
chr_multiplier            = scenario_tab['chr_multiplier', epi_scenario]
Itime_adav_chr_multiplier = scenario_tab['Itime_adav_chr_multiplier', epi_scenario]

# TOGGLE EFFECT OF AV ON TRANSMISSIBILITY
#if(av_ie_scenario == 0 ){av_transmissibility_adj = 1}
#if(av_ie_scenario == 50){av_transmissibility_adj = 0.5}

# TOGGLE BASELINE TRANSMISSIBILITY OF SEASON
#if(t_scenario == 'High'  ){sar_adj_noav = .086}
#if(t_scenario == 'Medium'){sar_adj_noav = .085}
#if(t_scenario == 'Low'   ){sar_adj_noav = .084}

# CALC TRANSMISSIBILITY INPUTS BASED ON SELECTIONS
sar_adj_os = sar_adj_noav*av_transmissibility_adj
SAR_fromage_noav = SAR_flutes*sar_adj_noav
SAR_fromage_os = SAR_flutes*sar_adj_os

#formula: beta = (num contacts per unit * prob transmission per contact) / N
betas_noav = (SAR_fromage_noav * contacts_combined) / age_pop_size #ops broadcast by column
betas_ad = (SAR_fromage_os * contacts_combined) / age_pop_size
betas_late = matrix(0, 5, 5) #betas_noav
betas_part = betas_noav


if(epi_scenario %in% c('lowtrans_baseline', 'hightrans_baseline')){
  betas_ad = betas_noav
}

##################################################
# CARE SEEKING #
##################################################

# from sineads appendix
p_seekcare = c(runif(1, 0.25, 0.35), 
               runif(1, 0.25, 0.35), 
               runif(1, 0.35, 0.50), 
               runif(1, 0.45, 0.55), 
               runif(1, 0.45, 0.55))     # non high risk
p_seekcare_hr = p_seekcare * 1.27        # high risk (1.27 is hr to std risk ratio from matt's paper, similar to 3 from sinead's appendix)

p_hr = c(0.05, 0.10, 0.20, 0.35, 0.55) # prob of high risk by age group
p_lr = 1-p_hr                          # prob not high risk by age group

p_seekcare_avg = p_seekcare * p_lr + p_seekcare_hr * p_hr  #avg probability of seeking care, including low and high risk 

p_nocare = 1-p_seekcare_avg

#p_ontime_given_seekcare = .34 # from matt's paper https://pmc.ncbi.nlm.nih.gov/articles/PMC4610008/table/T5/
p_ontime_given_seekcare = c(runif(1, 0.45, 0.55),
                            runif(1, 0.45, 0.55),
                            runif(1, 0.45, 0.55), 
                            runif(1, 0.40, 0.50),
                            runif(1, 0.35, 0.45))
p_late_given_seekcare = 1 - p_ontime_given_seekcare                       # about the same for lr versus hr

###############################################
# CARE SEEKING PROBABILITIES #
###############################################
p_adherence = .5

p_adav_given_ontime = .55 * p_adherence
p_part_given_ontime = .55 * (1- p_adherence)
p_noav_given_ontime = 1 - (p_adav_given_ontime + p_part_given_ontime)
  
p_av_given_late = .35
p_noav_given_late = 1 - p_av_given_late
  
p_time_adav = p_adav_given_ontime * p_ontime_given_seekcare * p_seekcare_avg
p_time_part = p_part_given_ontime * p_ontime_given_seekcare * p_seekcare_avg 
p_time_noav = p_noav_given_ontime * p_ontime_given_seekcare * p_seekcare_avg

p_late_av = p_av_given_late * p_late_given_seekcare * p_seekcare_avg
p_late_noav = p_noav_given_late * p_late_given_seekcare * p_seekcare_avg
  
p_nocare_noav = p_nocare


################################################
# RHOS - Isym to treatment classes #
################################################
pd_to_ontime_care = rnormTrunc(1, mean = 0.7 + 1.4 , sd = .2, min = 0, max = 2) #runif(1, 0.2, 2) #1.5
pd_to_late_care = rnorm(1, 2.7, .4) #runif(1, 2, 4)     #2.9
pd_to_no_care = 2 #runif(1, 0.2, 2.3)   #dummy

rho_time_adav = p_time_adav / pd_to_ontime_care
rho_time_part = p_time_part / pd_to_ontime_care
rho_time_noav = p_time_noav / pd_to_ontime_care

rho_late_av = p_late_av / pd_to_late_care
rho_late_noav = p_late_noav / pd_to_late_care

rho_nocare_noav = p_nocare_noav / pd_to_no_care

############################################### 
# RHO - Iasym to R #
###############################################
pd_asym_to_r = 2.7 #runif(1, 1.9, 3.7) #3
rho_asym_r = 1/pd_asym_to_r

###############################################
# MUs - treatment classes to recovered #
###############################################
pd_inf_post_ontime = .6 #max(0.01, 2.7 - pd_to_ontime_care) #runif(1, 1.7, 2.4)
pd_inf_post_late = 0.01 #runif(1, 0.01, 0.4) #0.1
pd_inf_post_nocare = 1 #runif(1, 1.7, 2.4) #1.5

mu_time_adav = 1/pd_inf_post_ontime  #these three mu's do not need to be same
mu_time_part = 1/pd_inf_post_ontime
mu_time_noav = 1/pd_inf_post_ontime

mu_late_av = 1/pd_inf_post_late
mu_late_noav = 1/pd_inf_post_late
mu_nocare_noav = 1/pd_inf_post_nocare


#############################################################################
# SET CASE HOSP RATIO AND CASE FATALITY RATIO #
#############################################################################

chr_by_age = c(143.44, 364.71, 178.16, 94.3, 11)                                #FROM FLU BURDEN DASH
deathhr_by_age = c(0.0132185, 0.0099303, 0.0294699, 0.0594536, 0.080231)        #FROM FLU BURDEN DASH

# TOGGLE BASELINE HOSP AND MORTALITY RATES
#if(s_scenario == 'Mild'  ){chr_multiplier = 1}
#if(s_scenario == 'Severe'){chr_multiplier = 1.5}

# CALC CHR AND CFR
chr_by_age_use = chr_multiplier*chr_by_age
deathhr_by_age_use = chr_multiplier*deathhr_by_age

#CHR MULTIPLIERS BY TREATMENT
#Itime_adav_chr_multiplier = 1                                                 #sinead varies from .2 to .75 (this is 1-that)
Itime_part_chr_multiplier = 1
Itime_noav_chr_multiplier = 1

Ilate_av_chr_multiplier = 1
Ilate_noav_chr_multiplier = 1

Inocare_noav_chr_multiplier = 1

#########################################
# set asymptomatic proportion
########################################
pd_latent = 1.5

p_asym = .3
p_sym = 1-p_asym

gamma_sym = p_sym / pd_latent 
gamma_asym = p_asym / pd_latent
