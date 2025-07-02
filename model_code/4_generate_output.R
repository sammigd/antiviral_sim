#DETERMINISTIC COMPARTMENTAL SEIR MODEL

#install.packages('deSolve')
library(deSolve)
library(tidyverse)
library(readxl)
library(here)
library(EnvStats)

set.seed(99)

#source(here('model_code/0_create_contact_matrix.R'))
source(here('model_code/1_set_fixed_params.R'))
source(here('model_code/3_create_model_function.R'))

scenario_tab = read.csv('inputs/scenarios.csv', row.names = 1)
epi_scenarios = names(scenario_tab)

#enumerate scenarios of interest
#severity_scenarios = c('Mild', 'Severe')
#antiviral_scenarios = c(0, # baseline - noav
                        #1)
#transmis_scenarios = c('Low', 'High')
#av_ie_scenarios = c(0, 50)

nsim = 10
#loop over scenarios of interest
res = data.frame(matrix(NA_real_, nrow = 0, ncol = 7))
res_ecurve = expand.grid(1:nsim, times, epi_scenarios)
names(res_ecurve) = c('iter', 'times', 'epi_scenario')

for (epi_scenario in epi_scenarios){
  #epi_scenario = epi_scenarios[1]
  for(i in 1:nsim){
    #i = 1   
    source('model_code/2_set_scenario_params.R')                     #set parameters matching the loop 
    tmp_params = c(betas_noav, betas_av, betas_three,     #generate input list for the model fn
                   gamma_sym, gamma_asym,
                   rho_asym_r, rho_r, 
                   rho_avzero, rho_avone, rho_avtwo, rho_avthree,
                   mu_avzero, mu_avone, mu_avtwo, mu_avthree)
  
    #run model
    tmp_out = as.data.frame(ode(y=state0, times=times, func=flu.model, parms=tmp_params, method="lsodes"))
    tmp_out_clean = tmp_out %>% 
      mutate(S_all = S1 + S2 + S3 + S4 + S5,
             E_all = E1 + E2 + E3 + E4 + E5,
             R_all = R1 + R2 + R3 + R4 + R5)
    
    res_ecurve$E_all[res_ecurve$iter == i &
                       res_ecurve$epi_scenario == epi_scenario] <- tmp_out_clean$E_all
    
    #calc final size
    tmp_final_size = max(tmp_out_clean$R_all) #the biggest / final recovered count
    
    tmp_out_clean$count_all_to_R = apply(tmp_out_clean[,c(47:76)], 1, sum)
    tmp_out_clean$count_sym_to_R = apply(tmp_out_clean[,c(47:71)], 1, sum) #these are already cumulative
    tmp_out_clean$count_asym_to_R = apply(tmp_out_clean[,c(72:76)], 1, sum)
    tmp_out_clean$n = apply(tmp_out_clean[,c(2:46)], 1, sum)
    
    sym_final_size = max(tmp_out_clean$count_sym_to_R)
    
    #calc hosp and mort counts
    for(age in 1:5){
      for (av in c('Iavzero', 'Iavone', 'Iavtwo', 'Iavthree', 'Isym_udx')){
        hosp_var_name = paste0('hosp_', av, age)               # colname for new hosp count columns
        death_var_name = paste0('death_', av, age)
        count_var_name = paste0('count_', av, '_to_R', age)    # colname for existing column with cumu transitions from I to R
        multiplier_use = paste0(av, '_chr_multiplier')         # multiplier for chr tagged to av treatment group
        chr_use = chr_by_age_use[age]
        
        tmp_out_clean[hosp_var_name] = (tmp_out_clean[count_var_name] / chr_use) * eval(parse(text = multiplier_use))
      }
    }
    
    total_hosp_a <- tmp_out_clean %>% summarise(across(hosp_Iavzero1:hosp_Isym_udx5, max))
    total_hosp_b <- sum(total_hosp_a) #sum across ages
  
    #deaths
    total_death = total_hosp_a
    for(age in 1:5){
      total_death[,str_detect(names(total_death), as.character(age))] <- total_death[,str_detect(names(total_death), as.character(age))] * deathhr_by_age_use[age]
    }
    total_death_b = sum(total_death) #sum across ages
    
    # total av - no longer accurate bc we arent counting ppl who seek care after day 3 who got av
    total_av = sum(tmp_out_clean %>% summarise(across(count_Iavzero_to_R1:count_Iavthree_to_R5, max)))
    #populate df w results
    tmp_results = c(i, epi_scenario, tmp_final_size, sym_final_size, total_hosp_b, total_death_b, total_av)
    res = rbind(res, tmp_results)
  }
}

sim_summary <- function(x){
  m = prettyNum(round(mean(x), 0), big.mark = ',')
  lb = prettyNum(round(quantile(x, 0.025),0), big.mark = ',')
  ub = prettyNum(round(quantile(x, 0.975),0), big.mark = ',')
  clean = paste0(m, ' (', lb, ' - ', ub, ')')
}

names(res) <- c('Iter', 'EpiScenario', 'FinalSize', 'SymFinalSize', 'HospCount', 'DeathCount', 'AVCount')
clean_res <- res %>% 
  mutate(across(FinalSize:AVCount, as.numeric))%>% 
  group_by(EpiScenario) %>%
  summarise(across(c(FinalSize, SymFinalSize, HospCount, DeathCount, AVCount), sim_summary))
clean_res
write.csv(clean_res, 'outputs/resultstab.csv')


#library(gt)
#gtsave(clean_res %>% 
#         gt() %>%
#         opt_horizontal_padding(scale = 3), 
#       filename = 'outputs/resultstab.png')



###########################################################
# PLOT EPI CURVES
###########################################################
res_ecurve |>
  #mutate(AVScenario = case_when(
  #  AVScenario == 0 ~ 'Baseline - No AV',
  #  AVScenario == 1 ~ 'Os like 2010',
  #  AVScenario == 2 ~ 'Bx like 2010')) |>
  group_by(times, epi_scenario) %>%
  summarise(E_all = mean(E_all)) %>%
  mutate(transmission = str_extract(epi_scenario, "^[^_]+"),
         antiviral = str_extract(epi_scenario, "(?<=_).*")) %>%
  ggplot(aes(x = times, y = E_all, group = antiviral, colour = antiviral)) + 
  geom_line() + 
  facet_grid(cols = vars(transmission)) + 
  ylab('Exposed Count') + xlab('Day (of 365)') + 
  theme_bw() 
ggsave('outputs/ecurves.png', width = 8, height = 4)



