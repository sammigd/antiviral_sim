#DETERMINISTIC COMPARTMENTAL SEIR MODEL

#install.packages('deSolve')
library(deSolve)
library(tidyverse)
library(readxl)
library(here)

set.seed(99)

#source(here('model_code/0_create_contact_matrix.R'))
source(here('model_code/1_set_fixed_params.R'))
source(here('model_code/3_create_model_function.R'))


#enumerate scenarios of interest
severity_scenarios = c('Mild', 'Severe')
antiviral_scenarios = c(0, # baseline - noav
                        1)
transmis_scenarios = c('Low', 'High')
av_ie_scenarios = c(0, 50)

nsim = 10
#loop over scenarios of interest
res = data.frame(matrix(NA_real_, nrow = 0, ncol = 7))
res_ecurve = expand.grid(1:nsim, times, severity_scenarios, transmis_scenarios, antiviral_scenarios, av_ie_scenarios)
names(res_ecurve) = c('iter', 'times', 'severity', 'transmission', 'AVScenario', 'AVIEScenario')

for(i in 1:nsim){
  for(s_scenario in severity_scenarios){
    for(antiviral_scenario in antiviral_scenarios){
      for (t_scenario in transmis_scenarios){
        for(av_ie_scenario in av_ie_scenarios){
          
          if(F){
            s_scenario = 'Mild'
            antiviral_scenario = 0
            t_scenario = 'Low'
            av_ie_scenario = 0
          }
          
          source('model_code/2_set_scenario_params.R')                     #set parameters matching the loop 
          tmp_params = c(betas_noav, betas_late, betas_part, betas_ad,     #generate input list for the model fn
                         gamma_sym, gamma_asym,
                         rho_asym_r, 
                         rho_time_adav, rho_time_part, rho_time_noav,
                         rho_late_av, rho_late_noav, rho_nocare_noav,
                         mu_time_adav, mu_time_part, mu_time_noav,
                         mu_late_av, mu_late_noav, mu_nocare_noav)
      
          #run model
          tmp_out = as.data.frame(ode(y=state0, times=times, func=flu.model, parms=tmp_params, method="lsodes"))
          tmp_out_clean = tmp_out %>% 
            mutate(S_all = S1 + S2 + S3 + S4 + S5,
                   E_all = E1 + E2 + E3 + E4 + E5,
                   R_all = R1 + R2 + R3 + R4 + R5)
          
          res_ecurve$E_all[res_ecurve$iter == i &
                             res_ecurve$severity == s_scenario & 
                             res_ecurve$transmission == t_scenario &
                             res_ecurve$AVScenario == antiviral_scenario &
                             res_ecurve$AVIEScenario == av_ie_scenario] <- tmp_out_clean$E_all
          
          #calc final size
          tmp_final_size = max(tmp_out_clean$R_all) #the biggest / final recovered count
          
          tmp_out_clean$count_all_to_R = apply(tmp_out_clean[,c(57:91)], 1, sum)
          tmp_out_clean$count_sym_to_R = apply(tmp_out_clean[,c(57:86)], 1, sum) #these are already cumulative
          tmp_out_clean$count_asym_to_R = apply(tmp_out_clean[,c(87:91)], 1, sum)
          tmp_out_clean$n = apply(tmp_out_clean[,c(2:56)], 1, sum)
          
          sym_final_size = max(tmp_out_clean$count_sym_to_R)
          
          #calc hosp and mort counts
          for(age in 1:5){
            for (av in c('Itime_adav', 'Itime_part', 'Itime_noav', 'Ilate_av', 'Ilate_noav', 'Inocare_noav')){
              hosp_var_name = paste0('hosp_', av, age)               # colname for new hosp count columns
              death_var_name = paste0('death_', av, age)
              count_var_name = paste0('count_', av, age, '_to_R')    # colname for existing column with cumu transitions from I to R
              multiplier_use = paste0(av, '_chr_multiplier')         # multiplier for chr tagged to av treatment group
              chr_use = chr_by_age_use[age]
              
              tmp_out_clean[hosp_var_name] = (tmp_out_clean[count_var_name] / chr_use) * eval(parse(text = multiplier_use))
            }
          }
          
          total_hosp_a <- tmp_out_clean %>% summarise(across(hosp_Itime_adav1:hosp_Inocare_noav5, max))
          total_hosp_b <- sum(total_hosp_a) #sum across ages
      
          #deaths
          total_death = total_hosp_a
          for(age in 1:5){
            total_death[,str_detect(names(total_death), as.character(age))] <- total_death[,str_detect(names(total_death), as.character(age))] * deathhr_by_age_use[age]
          }
          total_death_b = sum(total_death) #sum across ages
          
          # total av
          total_av = sum(tmp_out_clean %>% summarise(across(count_Itime_adav1_to_R:count_Ilate_av5_to_R, max))) - 
                     sum(tmp_out_clean %>% summarise(across(count_Itime_noav1_to_R:count_Itime_noav5_to_R, max)))
          
          #populate df w results
          tmp_results = c(i, s_scenario, t_scenario, antiviral_scenario, av_ie_scenario, tmp_final_size, sym_final_size, total_hosp_b, total_death_b, total_av)
          res = rbind(res, tmp_results)
        }
      }
    }
  }
}


sim_summary <- function(x){
  m = prettyNum(round(mean(x), 0), big.mark = ',')
  lb = prettyNum(round(quantile(x, 0.025),0), big.mark = ',')
  ub = prettyNum(round(quantile(x, 0.975),0), big.mark = ',')
  clean = paste0(m, ' (', lb, ' - ', ub, ')')
}

names(res) <- c('Iter', 'Severity\nScenario', 'Transmissibility\nScenario', 'AVScenario', 'AVIEScenario', 'FinalSize', 'SymFinalSize', 'HospCount', 'DeathCount', 'AVCount')
clean_res <- res %>% 
  mutate(across(FinalSize:AVCount, as.numeric),
         #across(FinalSize:AVCount, ~round(., 0)),
         #across(FinalSize:AVCount, ~prettyNum(., big.mark = ',')),
         AVScenario = case_when(
           AVScenario == 0 ~ 'Baseline - No AV',
           AVScenario == 1 ~ 'Os like 2010'
         )) %>% 
  group_by(`Severity\nScenario`, `Transmissibility\nScenario`, AVScenario, AVIEScenario) %>%
  summarise(across(c(FinalSize, SymFinalSize, HospCount, DeathCount, AVCount), sim_summary))
clean_res
#library(gt)
gtsave(clean_res %>% 
         gt() %>%
         opt_horizontal_padding(scale = 3), 
       filename = 'outputs/resultstab.png')



###########################################################
# PLOT EPI CURVES
###########################################################
res_ecurve |>
  mutate(AVScenario = case_when(
    AVScenario == 0 ~ 'Baseline - No AV',
    AVScenario == 1 ~ 'Os like 2010',
    AVScenario == 2 ~ 'Bx like 2010')) |>
  ggplot(aes(x = times, y = E_all, group = AVScenario, colour = AVScenario)) + 
  geom_line() + 
  facet_grid(cols = vars(severity, AVIEScenario), rows = vars(transmission)) + 
  ylab('Exposed Count') + xlab('Day (of 365)') + 
  theme_bw() + 
  labs(caption = "Severity scenarios (mild / severe) don't affect overall case counts, just hospitalizations and deaths")
ggsave('outputs/ecurves.png', width = 8, height = 4)



