#DETERMINISTIC COMPARTMENTAL SEIR MODEL

#install.packages('deSolve')
library(deSolve)
library(tidyverse)
library(readxl)
library(here)

source(here('model_code/0_create_contact_matrix.R'))
source(here('model_code/1_set_fixed_params.R'))
source(here('model_code/3_create_model_function.R'))


#enumerate scenarios of interest
severity_scenarios = c('Mild', 'Severe')
antiviral_scenarios = c(0, # baseline - noav
                        1, # baseline - observed (os)
                        2) # baseline - baloxavir
transmis_scenarios = c('Low', 'High')

#loop over scenarios of interest
res = data.frame(matrix(NA_real_, nrow = 0, ncol = 7))
res_ecurve = expand.grid(times, severity_scenarios, transmis_scenarios, antiviral_scenarios)
names(res_ecurve) = c('times', 'severity', 'transmission', 'AVScenario')

for(s_scenario in severity_scenarios){
  for(antiviral_scenario in antiviral_scenarios){
    for (t_scenario in transmis_scenarios){
      source('model_code/2_set_scenario_params.R')                          #set parameters matching the loop 
      tmp_params = c(betas_noav, betas_late, betas_nonad, betas_ontime,     #generate input list for the model fn
                     gamma, 
                     rho_r, rho_noav, rho_late, rho_ontime, rho_nonad,
                     mu_noav1, mu_late1, mu_nonad1, mu_ontime1)
  
      #run model
      tmp_out = as.data.frame(ode(y=state0, times=times, func=flu.model, parms=tmp_params, method="lsodes", scenario = 'Mild'))
      tmp_out_clean = tmp_out %>% 
        mutate(S_all = S1 + S2 + S3 + S4 + S5,
               E_all = E1 + E2 + E3 + E4 + E5,
               R_all = R1 + R2 + R3 + R4 + R5)
      
      res_ecurve$E_all[res_ecurve$severity == s_scenario & 
                         res_ecurve$transmission == t_scenario &
                         res_ecurve$AVScenario == antiviral_scenario] <- tmp_out_clean$E_all
      
      #calc final size
      tmp_final_size = max(tmp_out_clean$R_all)
      
      #calc hosp and mort counts
      for(age in 1:5){
        for (av in c('Inoav', 'Ilate', 'Inonad', 'Iontime', 'Iudx')){
          hosp_var_name = paste0('hosp_', av, age)               # colname for new hosp count columns
          death_var_name = paste0('death_', av, age)
          count_var_name = paste0('count_', av,age, '_to_R')     # colname for existing column with cumu transitions from I to R
          multiplier_use = paste0(av, '_chr_multiplier')         # multiplier for chr tagged to av treatment group
          chr_use = chr_by_age_use[age]
          
          tmp_out_clean[hosp_var_name] = (tmp_out_clean[count_var_name] / chr_use) * eval(parse(text = multiplier_use))
        }
      }
      
      total_hosp_a <- tmp_out_clean %>% summarise(across(hosp_Inoav1:hosp_Iudx5, max))
      total_hosp_b <- sum(total_hosp_a) #sum across ages
  
      #deaths
      total_death = total_hosp_a
      for(age in 1:5){
        total_death[,str_detect(names(total_death), as.character(age))] <- total_death[,str_detect(names(total_death), as.character(age))] * deathhr_by_age_use[age]
      }
      total_death_b = sum(total_death) #sum across ages
      
      # total av
      total_av = sum(tmp_out_clean %>% summarise(across(count_Ilate1_to_R:count_Iontime5_to_R, max)))
      
      #populate df w results
      tmp_results = c(s_scenario, t_scenario, antiviral_scenario, tmp_final_size, total_hosp_b, total_death_b, total_av)
      res = rbind(res, tmp_results)
    }
  }
}

names(res) <- c('Severity\nScenario', 'Transmissibility\nScenario', 'AVScenario', 'FinalSize', 'HospCount', 'DeathCount', 'AVCount')
clean_res <- res %>% 
  mutate(across(FinalSize:AVCount, as.numeric),
         across(FinalSize:AVCount, ~round(., 0)),
         across(FinalSize:AVCount, ~prettyNum(., big.mark = ',')),
         AVScenario = case_when(
           AVScenario == 0 ~ 'Baseline - No AV',
           AVScenario == 1 ~ 'Os like 2010',
           AVScenario == 2 ~ 'Bx like 2010'
         ))
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
  facet_grid(cols = vars(severity), rows = vars(transmission)) + 
  ylab('Exposed Count') + xlab('Day (of 365)') + 
  theme_bw() + 
  labs(caption = "Severity scenarios (mild / severe) don't affect overall case counts, just hospitalizations and deaths")
ggsave('outputs/ecurves.png', width = 8, height = 4)



