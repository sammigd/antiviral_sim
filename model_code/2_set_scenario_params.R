

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

av_transmissibility_adj   = scenario_tab['av_transmissibility_adj', epi_scenario]   # effect of av on transmissibility
sar_adj_noav              = scenario_tab['sar_adj_noav', epi_scenario]              # toggles baseline: low v high trans season              
p_adherence               = scenario_tab['p_adherence', epi_scenario]               # adherence to av
chr_multiplier            = scenario_tab['chr_multiplier', epi_scenario]            # toggles baseline: mild versus severe flu season
Itime_adav_chr_multiplier = scenario_tab['Itime_adav_chr_multiplier', epi_scenario] # effect of av on CHR

# CALC TRANSMISSIBILITY INPUTS BASED ON SELECTIONS
sar_adj_os = sar_adj_noav*av_transmissibility_adj
SAR_fromage_noav = SAR_flutes*sar_adj_noav
SAR_fromage_os = SAR_flutes*sar_adj_os

#formula: beta = (num contacts per unit * prob transmission per contact) / N
betas_noav = (SAR_fromage_noav * contacts_combined) / age_pop_size #ops broadcast by column
betas_av = (SAR_fromage_os * contacts_combined) / age_pop_size
betas_three = betas_noav #for now, trt on day three has no effect. can change. 

if(epi_scenario %in% c('lowtrans_baseline', 'hightrans_baseline')){
  betas_ad = betas_noav
}

#########################################
# set asymptomatic proportion
########################################
pd_latent = 1.5

p_asym = .3
p_sym = 1-p_asym

gamma_sym = p_sym / pd_latent 
gamma_asym = p_asym / pd_latent

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

##############################################
# DAY OF CARE SEEKING PROBS #
##############################################
p_d0 = 0.035 # from USFluVE 2023-2024
p_d1 = 0.20
p_d2 = 0.21
p_d3 = 0.22
p_over3 = 1 - p_d0 - p_d1 - p_d2  - p_d3

p_adherence = 0.5
p_prescribed = 0.55

p_avzero  = p_seekcare_avg * p_d0 * p_prescribed * p_adherence
p_avone   = p_seekcare_avg * p_d1 * p_prescribed * p_adherence
p_avtwo   = p_seekcare_avg * p_d2 * p_prescribed * p_adherence
p_avthree = p_seekcare_avg * p_d3 * p_prescribed * p_adherence

p_sym_nocare     = (1 - p_seekcare_avg)
p_sym_latecare   = p_seekcare_avg * p_over3
p_sym_care_noav  = p_seekcare_avg * (p_d0 + p_d1 + p_d2 + p_d3) * (1-p_prescribed)
p_sym_care_nonad = p_seekcare_avg * (p_d0 + p_d1 + p_d2 + p_d3) * p_prescribed * (1-p_adherence)

#line below should add up to 1
#p_avzero + p_avone + p_avtwo + p_avthree + p_sym_nocare + p_sym_latecare + p_sym_care_noav + p_sym_care_nonad

################################################
# RHOS - Isym to treatment classes and R #
################################################
dur_pre_sym_inf = 0.7
dur_inf = 3 + dur_pre_sym_inf
pd_Isym_to_avzero = 0.5 + dur_pre_sym_inf
pd_Isym_to_avone = 1 + dur_pre_sym_inf
pd_Isym_to_avtwo = 2 + dur_pre_sym_inf
pd_Isym_to_avthree = 2.5

rho_avzero  = p_avzero  / pd_Isym_to_avzero
rho_avone   = p_avone   / pd_Isym_to_avone
rho_avtwo   = p_avtwo   / pd_Isym_to_avtwo
rho_avthree = p_avthree / pd_Isym_to_avthree

rho_r = (p_sym_nocare + p_sym_latecare + p_sym_care_noav + p_sym_care_nonad) / 3

############################################### 
# RHO - Iasym #
###############################################
pd_asym_to_r = dur_inf 
rho_asym_r = 1/pd_asym_to_r # the same for each age group

###############################################
# MUs - treatment classes to R #
###############################################
mu_avzero  = rep(1 / (dur_inf - pd_Isym_to_avzero), 5) #the same for each age group
mu_avone   = rep(1 / (dur_inf - pd_Isym_to_avone), 5)
mu_avtwo   = rep(1 / (dur_inf - pd_Isym_to_avtwo), 5)
mu_avthree = rep(1 / (dur_inf - pd_Isym_to_avthree), 5)


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
Iavzero_chr_multiplier = 1
Iavone_chr_multiplier = 1

Iavtwo_chr_multiplier = 1
Iavthree_chr_multiplier = 1

Isym_udx_chr_multiplier = 1


