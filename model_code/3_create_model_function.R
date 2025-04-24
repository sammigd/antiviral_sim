##################################################
# WRITE MODEL FUNCTION ###########################
##################################################

flu.model <- function(t, state, parameters){
  with(as.list(c(state,parameters)),{
    
    #add seasonality to betas
    betas_noav = betas_noav * (1 + .1*sin(2*pi*(t/365)))
    betas_ad   = betas_ad   * (1 + .1*sin(2*pi*(t/365)))
    betas_late = betas_late * (1 + .1*sin(2*pi*(t/365)))
    betas_part = betas_part * (1 + .1*sin(2*pi*(t/365)))
    
    #TOTAL POPULATION - SPLIT INTO 5 AGE GROUPS
    #             |     PREDIAGNOSIS      |        TIMELY CARE SEEKERS              |     LATE CARE SEEKERS   | NO CARE SOUGHT| RECOVERED  
    N <- S1 + E1 + Isym_udx1 + Iasym_udx1 + Itime_adav1 + Itime_part1 + Itime_noav1 + Ilate_av1 + Ilate_noav1 + Inocare_noav1 + R1 +
         S2 + E2 + Isym_udx2 + Iasym_udx2 + Itime_adav2 + Itime_part2 + Itime_noav2 + Ilate_av2 + Ilate_noav2 + Inocare_noav2 + R2 +
         S3 + E3 + Isym_udx3 + Iasym_udx3 + Itime_adav3 + Itime_part3 + Itime_noav3 + Ilate_av3 + Ilate_noav3 + Inocare_noav3 + R3 +
         S4 + E4 + Isym_udx4 + Iasym_udx4 + Itime_adav4 + Itime_part4 + Itime_noav4 + Ilate_av4 + Ilate_noav4 + Inocare_noav4 + R4 +
         S5 + E5 + Isym_udx5 + Iasym_udx5 + Itime_adav5 + Itime_part5 + Itime_noav5 + Ilate_av5 + Ilate_noav5 + Inocare_noav5 + R5

    vec_Iall_noav = c(Iasym_udx1 + Isym_udx1 + Itime_noav1 + Ilate_noav1 + Inocare_noav1,
                      Iasym_udx2 + Isym_udx2 + Itime_noav2 + Ilate_noav2 + Inocare_noav2,
                      Iasym_udx3 + Isym_udx3 + Itime_noav3 + Ilate_noav3 + Inocare_noav3,
                      Iasym_udx4 + Isym_udx4 + Itime_noav4 + Ilate_noav4 + Inocare_noav4,
                      Iasym_udx5 + Isym_udx5 + Itime_noav5 + Ilate_noav5 + Inocare_noav5)
    vec_Itime_adav = c(Itime_adav1, Itime_adav2, Itime_adav3, Itime_adav4, Itime_adav5)
    vec_Itime_partav = c(Itime_part1, Itime_part2, Itime_part3, Itime_part4, Itime_part5) 
    vec_Ilate_av = c(Ilate_av1, Ilate_av2, Ilate_av3, Ilate_av4, Ilate_av5)
    
    dS1 <- -sum(betas_noav[1,]*S1*vec_Iall_noav) - sum(betas_ad[1,]*S1*vec_Itime_adav) - sum(betas_late[1,]*S1*vec_Ilate_av) - sum(betas_part[1,]*S1*vec_Itime_partav)
    dS2 <- -sum(betas_noav[2,]*S2*vec_Iall_noav) - sum(betas_ad[2,]*S2*vec_Itime_adav) - sum(betas_late[2,]*S2*vec_Ilate_av) - sum(betas_part[2,]*S2*vec_Itime_partav)
    dS3 <- -sum(betas_noav[3,]*S3*vec_Iall_noav) - sum(betas_ad[3,]*S3*vec_Itime_adav) - sum(betas_late[3,]*S3*vec_Ilate_av) - sum(betas_part[3,]*S3*vec_Itime_partav)
    dS4 <- -sum(betas_noav[4,]*S4*vec_Iall_noav) - sum(betas_ad[4,]*S4*vec_Itime_adav) - sum(betas_late[4,]*S4*vec_Ilate_av) - sum(betas_part[4,]*S4*vec_Itime_partav)
    dS5 <- -sum(betas_noav[5,]*S5*vec_Iall_noav) - sum(betas_ad[5,]*S5*vec_Itime_adav) - sum(betas_late[5,]*S5*vec_Ilate_av) - sum(betas_part[5,]*S5*vec_Itime_partav)
    
    dE1 <-  sum(betas_noav[1,]*S1*vec_Iall_noav) + sum(betas_ad[1,]*S1*vec_Itime_adav) + sum(betas_late[1,]*S1*vec_Ilate_av) + sum(betas_part[1,]*S1*vec_Itime_partav) - gamma_asym*E1 - gamma_sym*E1
    dE2 <-  sum(betas_noav[2,]*S2*vec_Iall_noav) + sum(betas_ad[2,]*S2*vec_Itime_adav) + sum(betas_late[2,]*S2*vec_Ilate_av) + sum(betas_part[2,]*S2*vec_Itime_partav) - gamma_asym*E2 - gamma_sym*E2
    dE3 <-  sum(betas_noav[3,]*S3*vec_Iall_noav) + sum(betas_ad[3,]*S3*vec_Itime_adav) + sum(betas_late[3,]*S3*vec_Ilate_av) + sum(betas_part[3,]*S3*vec_Itime_partav) - gamma_asym*E3 - gamma_sym*E3
    dE4 <-  sum(betas_noav[4,]*S4*vec_Iall_noav) + sum(betas_ad[4,]*S4*vec_Itime_adav) + sum(betas_late[4,]*S4*vec_Ilate_av) + sum(betas_part[4,]*S4*vec_Itime_partav) - gamma_asym*E4 - gamma_sym*E4
    dE5 <-  sum(betas_noav[5,]*S5*vec_Iall_noav) + sum(betas_ad[5,]*S5*vec_Itime_adav) + sum(betas_late[5,]*S5*vec_Ilate_av) + sum(betas_part[5,]*S5*vec_Itime_partav) - gamma_asym*E5 - gamma_sym*E5
    
    dIasym_udx1 = gamma_asym*E1 - rho_asym_r*Iasym_udx1                         # undiagnosed asymptomatic compartments
    dIasym_udx2 = gamma_asym*E2 - rho_asym_r*Iasym_udx2 
    dIasym_udx3 = gamma_asym*E3 - rho_asym_r*Iasym_udx3 
    dIasym_udx4 = gamma_asym*E4 - rho_asym_r*Iasym_udx4 
    dIasym_udx5 = gamma_asym*E5 - rho_asym_r*Iasym_udx5 
    
    dIsym_udx1 <- gamma_sym*E1 - Isym_udx1*(rho_time_adav[1] + rho_time_part[1] + rho_time_noav[1] + rho_late_av[1] + rho_late_noav[1] + rho_nocare_noav[1])      # udx symptomatic compartments
    dIsym_udx2 <- gamma_sym*E2 - Isym_udx2*(rho_time_adav[2] + rho_time_part[2] + rho_time_noav[2] + rho_late_av[2] + rho_late_noav[2] + rho_nocare_noav[2])      # from: E
    dIsym_udx3 <- gamma_sym*E3 - Isym_udx3*(rho_time_adav[3] + rho_time_part[3] + rho_time_noav[3] + rho_late_av[3] + rho_late_noav[3] + rho_nocare_noav[3])      # to: timely care: adherent, partial, noav
    dIsym_udx4 <- gamma_sym*E4 - Isym_udx4*(rho_time_adav[4] + rho_time_part[4] + rho_time_noav[4] + rho_late_av[4] + rho_late_noav[4] + rho_nocare_noav[4])      #     late care: av, noav
    dIsym_udx5 <- gamma_sym*E5 - Isym_udx5*(rho_time_adav[5] + rho_time_part[5] + rho_time_noav[5] + rho_late_av[5] + rho_late_noav[5] + rho_nocare_noav[5])      #     no care: noav
    
    dItime_adav1 = rho_time_adav[1]*Isym_udx1 - mu_time_adav*Itime_adav1     # timely careseeking - adherent antivirals
    dItime_adav2 = rho_time_adav[2]*Isym_udx2 - mu_time_adav*Itime_adav2
    dItime_adav3 = rho_time_adav[3]*Isym_udx3 - mu_time_adav*Itime_adav3
    dItime_adav4 = rho_time_adav[4]*Isym_udx4 - mu_time_adav*Itime_adav4
    dItime_adav5 = rho_time_adav[5]*Isym_udx5 - mu_time_adav*Itime_adav5
    
    dItime_part1 = rho_time_part[1]*Isym_udx1 - mu_time_part*Itime_part1     # timely careseeking - partially adherent antivirals
    dItime_part2 = rho_time_part[2]*Isym_udx2 - mu_time_part*Itime_part2
    dItime_part3 = rho_time_part[3]*Isym_udx3 - mu_time_part*Itime_part3
    dItime_part4 = rho_time_part[4]*Isym_udx4 - mu_time_part*Itime_part4
    dItime_part5 = rho_time_part[5]*Isym_udx5 - mu_time_part*Itime_part5
    
    dItime_noav1 = rho_time_noav[1]*Isym_udx1 - mu_time_noav*Itime_noav1     # timely careseeking - no antivirals
    dItime_noav2 = rho_time_noav[2]*Isym_udx2 - mu_time_noav*Itime_noav2
    dItime_noav3 = rho_time_noav[3]*Isym_udx3 - mu_time_noav*Itime_noav3
    dItime_noav4 = rho_time_noav[4]*Isym_udx4 - mu_time_noav*Itime_noav4
    dItime_noav5 = rho_time_noav[5]*Isym_udx5 - mu_time_noav*Itime_noav5
    
    dIlate_av1 = rho_late_av[1]*Isym_udx1 - mu_late_av*Ilate_av1             # late careseeking - antivirals
    dIlate_av2 = rho_late_av[2]*Isym_udx2 - mu_late_av*Ilate_av2
    dIlate_av3 = rho_late_av[3]*Isym_udx3 - mu_late_av*Ilate_av3
    dIlate_av4 = rho_late_av[4]*Isym_udx4 - mu_late_av*Ilate_av4
    dIlate_av5 = rho_late_av[5]*Isym_udx5 - mu_late_av*Ilate_av5
    
    dIlate_noav1 = rho_late_noav[1]*Isym_udx1 - mu_late_noav*Ilate_noav1     # late careseeking - no antivirals
    dIlate_noav2 = rho_late_noav[2]*Isym_udx2 - mu_late_noav*Ilate_noav2
    dIlate_noav3 = rho_late_noav[3]*Isym_udx3 - mu_late_noav*Ilate_noav3
    dIlate_noav4 = rho_late_noav[4]*Isym_udx4 - mu_late_noav*Ilate_noav4
    dIlate_noav5 = rho_late_noav[5]*Isym_udx5 - mu_late_noav*Ilate_noav5
    
    dInocare_noav1 = rho_nocare_noav[1]*Isym_udx1 - mu_nocare_noav*Inocare_noav1     # no careseeking - no antivirals
    dInocare_noav2 = rho_nocare_noav[2]*Isym_udx2 - mu_nocare_noav*Inocare_noav2
    dInocare_noav3 = rho_nocare_noav[3]*Isym_udx3 - mu_nocare_noav*Inocare_noav3
    dInocare_noav4 = rho_nocare_noav[4]*Isym_udx4 - mu_nocare_noav*Inocare_noav4
    dInocare_noav5 = rho_nocare_noav[5]*Isym_udx5 - mu_nocare_noav*Inocare_noav5
    
    dR1 <- mu_time_adav*Itime_adav1 + mu_time_part*Itime_part1 +   mu_time_noav*Itime_noav1 + 
             mu_late_av*Ilate_av1   + mu_late_noav*Ilate_noav1 + mu_nocare_noav*Inocare_noav1 + rho_asym_r*Iasym_udx1
    dR2 <- mu_time_adav*Itime_adav2 + mu_time_part*Itime_part2 +   mu_time_noav*Itime_noav2 + 
             mu_late_av*Ilate_av2   + mu_late_noav*Ilate_noav2 + mu_nocare_noav*Inocare_noav2 + rho_asym_r*Iasym_udx2
    dR3 <- mu_time_adav*Itime_adav3 + mu_time_part*Itime_part3 +   mu_time_noav*Itime_noav3 + 
             mu_late_av*Ilate_av3   + mu_late_noav*Ilate_noav3 + mu_nocare_noav*Inocare_noav3 + rho_asym_r*Iasym_udx3
    dR4 <- mu_time_adav*Itime_adav4 + mu_time_part*Itime_part4 +   mu_time_noav*Itime_noav4 + 
             mu_late_av*Ilate_av4   + mu_late_noav*Ilate_noav4 + mu_nocare_noav*Inocare_noav4 + rho_asym_r*Iasym_udx4
    dR5 <- mu_time_adav*Itime_adav5 + mu_time_part*Itime_part5 +   mu_time_noav*Itime_noav5 + 
             mu_late_av*Ilate_av5   + mu_late_noav*Ilate_noav5 + mu_nocare_noav*Inocare_noav5 + rho_asym_r*Iasym_udx5
    
    #COUNT TRANSITIONS FROM EACH COMPARTMENT TO RECOVERED
    count_Itime_adav1_to_R = mu_time_adav*Itime_adav1 #timely care seekers
    count_Itime_adav2_to_R = mu_time_adav*Itime_adav2
    count_Itime_adav3_to_R = mu_time_adav*Itime_adav3
    count_Itime_adav4_to_R = mu_time_adav*Itime_adav4
    count_Itime_adav5_to_R = mu_time_adav*Itime_adav5
    
    count_Itime_part1_to_R = mu_time_part*Itime_part1
    count_Itime_part2_to_R = mu_time_part*Itime_part2
    count_Itime_part3_to_R = mu_time_part*Itime_part3
    count_Itime_part4_to_R = mu_time_part*Itime_part4
    count_Itime_part5_to_R = mu_time_part*Itime_part5
    
    count_Itime_noav1_to_R = mu_time_noav*Itime_noav1
    count_Itime_noav2_to_R = mu_time_noav*Itime_noav2
    count_Itime_noav3_to_R = mu_time_noav*Itime_noav3
    count_Itime_noav4_to_R = mu_time_noav*Itime_noav4
    count_Itime_noav5_to_R = mu_time_noav*Itime_noav5
    
    count_Ilate_av1_to_R = mu_late_av*Ilate_av1        #late careseekers
    count_Ilate_av2_to_R = mu_late_av*Ilate_av2
    count_Ilate_av3_to_R = mu_late_av*Ilate_av3
    count_Ilate_av4_to_R = mu_late_av*Ilate_av4
    count_Ilate_av5_to_R = mu_late_av*Ilate_av5
    
    count_Ilate_noav1_to_R = mu_late_noav*Ilate_noav1
    count_Ilate_noav2_to_R = mu_late_noav*Ilate_noav2
    count_Ilate_noav3_to_R = mu_late_noav*Ilate_noav3
    count_Ilate_noav4_to_R = mu_late_noav*Ilate_noav4
    count_Ilate_noav5_to_R = mu_late_noav*Ilate_noav5
    
    count_Inocare_noav1_to_R = mu_nocare_noav*Inocare_noav1 # no care   
    count_Inocare_noav2_to_R = mu_nocare_noav*Inocare_noav2
    count_Inocare_noav3_to_R = mu_nocare_noav*Inocare_noav3
    count_Inocare_noav4_to_R = mu_nocare_noav*Inocare_noav4
    count_Inocare_noav5_to_R = mu_nocare_noav*Inocare_noav5
    
    count_Iasym_udx1_to_R = rho_asym_r*Iasym_udx1
    count_Iasym_udx2_to_R = rho_asym_r*Iasym_udx2
    count_Iasym_udx3_to_R = rho_asym_r*Iasym_udx3
    count_Iasym_udx4_to_R = rho_asym_r*Iasym_udx4
    count_Iasym_udx5_to_R = rho_asym_r*Iasym_udx5
    
    

    return(list(c(dS1, dS2, dS3, dS4, dS5, 
                  dE1, dE2, dE3, dE4, dE5,
                  dIsym_udx1, dIsym_udx2, dIsym_udx3, dIsym_udx4, dIsym_udx4,
                  dIasym_udx1, dIasym_udx2, dIasym_udx3, dIasym_udx4, dIasym_udx4,
                  dItime_adav1, dItime_adav2, dItime_adav3, dItime_adav4, dItime_adav5,
                  dItime_part1, dItime_part2, dItime_part3, dItime_part4, dItime_part5,
                  dItime_noav1, dItime_noav2, dItime_noav3, dItime_noav4, dItime_noav5,
                  dIlate_av1, dIlate_av2, dIlate_av3, dIlate_av4, dIlate_av5,
                  dIlate_noav1, dIlate_noav2, dIlate_noav3, dIlate_noav4, dIlate_noav5,
                  dInocare_noav1, dInocare_noav2, dInocare_noav3, dInocare_noav4, dInocare_noav5,
                  dR1, dR2, dR3, dR4, dR5,
                  count_Itime_adav1_to_R, count_Itime_adav2_to_R, count_Itime_adav3_to_R, count_Itime_adav4_to_R, count_Itime_adav5_to_R,
                  count_Itime_part1_to_R, count_Itime_part2_to_R, count_Itime_part3_to_R, count_Itime_part4_to_R, count_Itime_part5_to_R,
                  count_Itime_noav1_to_R, count_Itime_noav2_to_R, count_Itime_noav3_to_R, count_Itime_noav4_to_R, count_Itime_noav5_to_R,
                  count_Ilate_av1_to_R, count_Ilate_av2_to_R, count_Ilate_av3_to_R, count_Ilate_av4_to_R, count_Ilate_av5_to_R,
                  count_Ilate_noav1_to_R, count_Ilate_noav2_to_R, count_Ilate_noav3_to_R, count_Ilate_noav4_to_R, count_Ilate_noav5_to_R,
                  count_Inocare_noav1_to_R, count_Inocare_noav2_to_R, count_Inocare_noav3_to_R, count_Inocare_noav4_to_R, count_Inocare_noav5_to_R,
                  count_Iasym_udx1_to_R, count_Iasym_udx2_to_R, count_Iasym_udx3_to_R, count_Iasym_udx4_to_R, count_Iasym_udx5_to_R)))
  })
}
