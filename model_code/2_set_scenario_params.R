

##################################################
# CALC BETAS #
##################################################

#SAR HARDCODED IN FROM UNPUBLISHED SAR IMPUTATION WORK
#for age groups 0-4, 5-17, 18-49, 50-64,50+
if(antiviral_scenario == 0){
  SAR_flutes = c(0.21, 0.14, 0.19, 0.25, 0.25)
}


if(t_scenario == 'High'){
  sar_adj_noav = .0835
  sar_adj_os = .0435
}

if(t_scenario == 'Low'){
  sar_adj_noav = .082
  sar_adj_os = .042
}

SAR_fromage_noav = SAR_flutes*sar_adj_noav

if(antiviral_scenario %in% c(0,1)){
  SAR_fromage_os = SAR_flutes*sar_adj_os
}
if(antiviral_scenario == 2){ #bx
  SAR_fromage_os = SAR_flutes*sar_adj_os*.5
}

#formula: beta = (num contacts per unit * prob transmission per contact) / N
betas_noav = (SAR_fromage_noav * contacts_combined) / age_pop_size #ops broadcast by column
betas_ontime = (SAR_fromage_os * contacts_combined) / age_pop_size
betas_late = betas_noav
betas_nonad = betas_noav



##################################################
# CARE SEEKING #
##################################################

                                         # from sineads appendix
p_seekcare = c(0.3, 0.3, .425, 0.5, 0.5) # non high risk
p_seekcare_hr = p_seekcare * 1.27        # high risk (1.27 is hr to std risk ratio from matt's paper, similar to 3 from sinead's appendix)

p_hr = c(0.05, 0.10, 0.20, 0.35, 0.55) # prob of high risk by age group
p_lr = 1-p_hr                          # prob not high risk by age group

p_seekcare_avg = p_seekcare * p_lr + p_seekcare_hr * p_hr  #avg probability of seeking care, including low and high risk 

p_nocare = 1-p_seekcare_avg

p_ontime_given_seekcare = .34 # from matt's paper https://pmc.ncbi.nlm.nih.gov/articles/PMC4610008/table/T5/
p_late_given_seekcare = 1 - p_ontime_given_seekcare                       # about the same for lr versus hr



#################################################
# AV PRESCRIBING PARAMETERS #
#################################################

if(antiviral_scenario == 0){ #baseline - no antivirals
  p_seekcare_anyav = p_seekcare * 0
  
  p_seekcare_ontime = p_seekcare * 0
  p_seekcare_nonad = p_seekcare * 0
  p_seekcare_late = p_seekcare * 0
  p_seekcare_noav = p_seekcare * 1
}

if(antiviral_scenario == 1){ #baseline - observed antivirals
  
  #set amount of nonadherence
  p_adherence_given_ontime = .5
  p_nonadherence_given_ontime = 1-p_adherence_given_ontime
  
  #set av prescribing of care seeking
  p_seekcare_anyav = .5 #p av | seek care
  
  #calculate probabilities for each treatment bin
  p_seekcare_ontime = p_seekcare_avg * p_seekcare_anyav * p_ontime_given_seekcare * p_adherence_given_ontime      # * p adherence
  p_seekcare_nonad =  p_seekcare_avg * p_seekcare_anyav * p_ontime_given_seekcare * p_nonadherence_given_ontime   # * p nonadherence
  p_seekcare_late =   p_seekcare_avg * p_late_given_seekcare * p_seekcare_anyav
  p_seekcare_noav =   p_seekcare_avg * (1- p_seekcare_anyav) #0.297
}

if(antiviral_scenario == 2){ #baseline - baloxavir
  
  #set amount of nonadherence - with bx 100% adherence bc one pill
  p_adherence_given_ontime = 1
  p_nonadherence_given_ontime = 1-p_adherence_given_ontime
  
  #set av prescribing of care seeking
  p_seekcare_anyav = .5 #p av | seek care
  
  #calculate probabilities for each treatment bin
  p_seekcare_ontime = p_seekcare_avg * p_seekcare_anyav * p_ontime_given_seekcare * p_adherence_given_ontime      # * p adherence
  p_seekcare_nonad =  p_seekcare_avg * p_seekcare_anyav * p_ontime_given_seekcare * p_nonadherence_given_ontime   # * p nonadherence
  p_seekcare_late =   p_seekcare_avg * p_late_given_seekcare * p_seekcare_anyav
  p_seekcare_noav =   p_seekcare_avg * (1- p_seekcare_anyav) 
}

#calc parms for Iudx to Idx
rho_noav = p_seekcare_noav / udx_to_care_pd
rho_late = p_seekcare_late / udx_to_care_pd
rho_nonad = p_seekcare_nonad / udx_to_care_pd
rho_ontime = p_seekcare_ontime / udx_to_care_pd
rho_r = p_nocare / udx_to_r_pd #dummy


#rates from care seeking infected to recovered
mu_ontime1 = 1/av_rec_pd
mu_nonad1 = 1/noav_rec_pd
mu_late1 = 1/noav_rec_pd
mu_noav1 = 1/noav_rec_pd

#############################################################################
# SET CASE HOSP RATIO AND CASE FATALITY RATIO #
#############################################################################

chr_by_age = c(143.44, 364.71, 178.16, 94.3, 11)                              #FROM FLU BURDEN DASH
deathhr_by_age = c(0.0132185, 0.0099303, 0.0294699, 0.0594536, 0.080231)      #FROM FLU BURDEN DASH

if(s_scenario == 'Mild'){
  chr_by_age_use = 1*chr_by_age
  deathhr_by_age_use = 1*deathhr_by_age
}
if(s_scenario == 'Severe'){
  chr_by_age_use = 1.5*chr_by_age
  deathhr_by_age_use = 1.5*deathhr_by_age
}

#CHR MULTIPLIERS BY TREATMENT
if(antiviral_scenario %in% c(0,1)){ #no av, os scenarios
  Inoav_chr_multiplier = 1
  Inonad_chr_multiplier = 1
  Ilate_chr_multiplier = 1
  Iontime_chr_multiplier = .8                                                   #sinead varies from .2 to .75 (this is 1-that)
  Iudx_chr_multiplier = 1
}



if(antiviral_scenario == 2){ #bx scenario
  Inoav_chr_multiplier = 1
  Inonad_chr_multiplier = 1
  Ilate_chr_multiplier = 1
  Iontime_chr_multiplier = .75                                                   #sinead varies from .2 to .75 (this is 1-that)
  Iudx_chr_multiplier = 1
}
