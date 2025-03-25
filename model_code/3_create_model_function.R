##################################################
# WRITE MODEL FUNCTION ###########################
##################################################

flu.model <- function(t, state, parameters, scenario){
  with(as.list(c(state,parameters)),{
    
    #add seasonality to betas
    betas_noav = betas_noav* (1 + .1*sin(2*pi*(t/365)))
    betas_ontime = betas_ontime* (1 + .1*sin(2*pi*(t/365)))
    betas_late = betas_late* (1 + .1*sin(2*pi*(t/365)))
    betas_nonad = betas_nonad* (1 + .1*sin(2*pi*(t/365)))
    
    #TOTAL POPULATION - SPLIT INTO 5 AGE GROUPS
    #START WITH ONE AGE GROUP
    N <- S1 + E1 + Iudx1 + Iontime1 + Ilate1 + Inonad1 +  Inoav1 + R1 + 
      S2 + E2 + Iudx2 + Iontime2 + Ilate2 + Inonad2 +  Inoav2 + R2 +
      S3 + E3 + Iudx3 + Iontime3 + Ilate3 + Inonad3 +  Inoav3 + R3 +
      S4 + E4 + Iudx4 + Iontime4 + Ilate4 + Inonad4 +  Inoav4 + R4 +
      S5 + E5 + Iudx5 + Iontime5 + Ilate5 + Inonad5 +  Inoav5 + R5
    
    vec_Iudx_or_Inoav = c(Iudx1+Inoav1, Iudx2+Inoav2, Iudx3+Inoav3, Iudx4+Inoav4, Iudx5+Inoav5)
    vec_Iontime = c(Iontime1, Iontime2, Iontime3, Iontime4, Iontime5)
    vec_Inonad = c(Inonad1, Inonad2, Inonad3, Inonad4, Inonad5)
    vec_Ilate = c(Ilate1, Ilate2, Ilate3, Ilate4, Ilate5)
    
    dS1 <- -sum(betas_noav[1,]*S1*vec_Iudx_or_Inoav) - sum(betas_ontime[1,]*S1*vec_Iontime) - sum(betas_late[1,]*S1*vec_Ilate) - sum(betas_nonad[1,]*S1*vec_Inonad)
    dS2 <- -sum(betas_noav[2,]*S2*vec_Iudx_or_Inoav) - sum(betas_ontime[2,]*S2*vec_Iontime) - sum(betas_late[2,]*S2*vec_Ilate) - sum(betas_nonad[2,]*S2*vec_Inonad)
    dS3 <- -sum(betas_noav[3,]*S3*vec_Iudx_or_Inoav) - sum(betas_ontime[3,]*S3*vec_Iontime) - sum(betas_late[3,]*S3*vec_Ilate) - sum(betas_nonad[3,]*S3*vec_Inonad)
    dS4 <- -sum(betas_noav[4,]*S4*vec_Iudx_or_Inoav) - sum(betas_ontime[4,]*S4*vec_Iontime) - sum(betas_late[4,]*S4*vec_Ilate) - sum(betas_nonad[4,]*S4*vec_Inonad)
    dS5 <- -sum(betas_noav[5,]*S5*vec_Iudx_or_Inoav) - sum(betas_ontime[5,]*S5*vec_Iontime) - sum(betas_late[5,]*S5*vec_Ilate) - sum(betas_nonad[5,]*S5*vec_Inonad)
    
    dE1 <- sum(betas_noav[1,]*S1*vec_Iudx_or_Inoav) + sum(betas_ontime[1,]*S1*vec_Iontime) + sum(betas_late[1,]*S1*vec_Ilate) + sum(betas_nonad[1,]*S1*vec_Inonad) - gamma*E1
    dE2 <- sum(betas_noav[2,]*S2*vec_Iudx_or_Inoav) + sum(betas_ontime[2,]*S2*vec_Iontime) + sum(betas_late[2,]*S2*vec_Ilate) + sum(betas_nonad[2,]*S2*vec_Inonad) - gamma*E2
    dE3 <- sum(betas_noav[3,]*S3*vec_Iudx_or_Inoav) + sum(betas_ontime[3,]*S3*vec_Iontime) + sum(betas_late[3,]*S3*vec_Ilate) + sum(betas_nonad[3,]*S3*vec_Inonad) - gamma*E3
    dE4 <- sum(betas_noav[4,]*S4*vec_Iudx_or_Inoav) + sum(betas_ontime[4,]*S4*vec_Iontime) + sum(betas_late[4,]*S4*vec_Ilate) + sum(betas_nonad[4,]*S4*vec_Inonad) - gamma*E4
    dE5 <- sum(betas_noav[5,]*S5*vec_Iudx_or_Inoav) + sum(betas_ontime[5,]*S5*vec_Iontime) + sum(betas_late[5,]*S5*vec_Ilate) + sum(betas_nonad[5,]*S5*vec_Inonad) - gamma*E5
    
    #intentionally keeping rho and mu constant across age groups for now
    dIudx1 <- gamma*E1 - rho_ontime[1]*Iudx1 - rho_late[1]*Iudx1 - rho_noav[1]*Iudx1 - rho_nonad[1]*Iudx1 - rho_r[1]*Iudx1
    dIudx2 <- gamma*E2 - rho_ontime[2]*Iudx2 - rho_late[2]*Iudx2 - rho_noav[2]*Iudx2 - rho_nonad[2]*Iudx2 - rho_r[2]*Iudx2
    dIudx3 <- gamma*E3 - rho_ontime[3]*Iudx3 - rho_late[3]*Iudx3 - rho_noav[3]*Iudx3 - rho_nonad[3]*Iudx3 - rho_r[3]*Iudx3
    dIudx4 <- gamma*E4 - rho_ontime[4]*Iudx4 - rho_late[4]*Iudx4 - rho_noav[4]*Iudx4 - rho_nonad[4]*Iudx4 - rho_r[4]*Iudx4
    dIudx5 <- gamma*E5 - rho_ontime[5]*Iudx5 - rho_late[5]*Iudx5 - rho_noav[5]*Iudx5 - rho_nonad[5]*Iudx5 - rho_r[5]*Iudx5
    
    dIontime1 <-  rho_ontime[1]*Iudx1 - mu_ontime1*Iontime1
    dIontime2 <-  rho_ontime[2]*Iudx2 - mu_ontime1*Iontime2
    dIontime3 <-  rho_ontime[3]*Iudx3 - mu_ontime1*Iontime3
    dIontime4 <-  rho_ontime[4]*Iudx4 - mu_ontime1*Iontime4
    dIontime5 <-  rho_ontime[5]*Iudx5 - mu_ontime1*Iontime5
    
    dIlate1 <-  rho_late[1]*Iudx1 - mu_late1*Ilate1
    dIlate2 <-  rho_late[2]*Iudx2 - mu_late1*Ilate2
    dIlate3 <-  rho_late[3]*Iudx3 - mu_late1*Ilate3
    dIlate4 <-  rho_late[4]*Iudx4 - mu_late1*Ilate4
    dIlate5 <-  rho_late[5]*Iudx5 - mu_late1*Ilate5
    
    dInonad1 <- rho_nonad[1]*Iudx1 - mu_nonad1*Inonad1
    dInonad2 <- rho_nonad[2]*Iudx2 - mu_nonad1*Inonad2
    dInonad3 <- rho_nonad[3]*Iudx3 - mu_nonad1*Inonad3
    dInonad4 <- rho_nonad[4]*Iudx4 - mu_nonad1*Inonad4
    dInonad5 <- rho_nonad[5]*Iudx5 - mu_nonad1*Inonad5
    
    dInoav1 <- rho_noav[1]*Iudx1 - mu_noav1*Inoav1
    dInoav2 <- rho_noav[2]*Iudx2 - mu_noav1*Inoav2
    dInoav3 <- rho_noav[3]*Iudx3 - mu_noav1*Inoav3
    dInoav4 <- rho_noav[4]*Iudx4 - mu_noav1*Inoav4
    dInoav5 <- rho_noav[5]*Iudx5 - mu_noav1*Inoav5
    
    dR1 <- mu_ontime1*Iontime1 + mu_late1*Ilate1 + mu_nonad1*Inonad1 + mu_noav1*Inoav1 + rho_r[1]*Iudx1
    dR2 <- mu_ontime1*Iontime2 + mu_late1*Ilate2 + mu_nonad1*Inonad2 + mu_noav1*Inoav2 + rho_r[2]*Iudx2
    dR3 <- mu_ontime1*Iontime3 + mu_late1*Ilate3 + mu_nonad1*Inonad3 + mu_noav1*Inoav3 + rho_r[3]*Iudx3
    dR4 <- mu_ontime1*Iontime4 + mu_late1*Ilate4 + mu_nonad1*Inonad4 + mu_noav1*Inoav4 + rho_r[4]*Iudx4
    dR5 <- mu_ontime1*Iontime5 + mu_late1*Ilate5 + mu_nonad1*Inonad5 + mu_noav1*Inoav5 + rho_r[5]*Iudx5
    
    #COUNT TRANSITIONS FROM EACH COMPARTMENT TO RECOVERED
    count_Inoav1_to_R = mu_noav1*Inoav1
    count_Inoav2_to_R = mu_noav1*Inoav2
    count_Inoav3_to_R = mu_noav1*Inoav3
    count_Inoav4_to_R = mu_noav1*Inoav4
    count_Inoav5_to_R = mu_noav1*Inoav5
    
    count_Ilate1_to_R = mu_late1*Ilate1
    count_Ilate2_to_R = mu_late1*Ilate2
    count_Ilate3_to_R = mu_late1*Ilate3
    count_Ilate4_to_R = mu_late1*Ilate4
    count_Ilate5_to_R = mu_late1*Ilate5
    
    count_Inonad1_to_R = mu_nonad1*Inonad1
    count_Inonad2_to_R = mu_nonad1*Inonad2
    count_Inonad3_to_R = mu_nonad1*Inonad3
    count_Inonad4_to_R = mu_nonad1*Inonad4
    count_Inonad5_to_R = mu_nonad1*Inonad5
    
    count_Iontime1_to_R = mu_ontime1*Iontime1
    count_Iontime2_to_R = mu_ontime1*Iontime2
    count_Iontime3_to_R = mu_ontime1*Iontime3
    count_Iontime4_to_R = mu_ontime1*Iontime4
    count_Iontime5_to_R = mu_ontime1*Iontime5
    
    count_Iudx1_to_R = rho_r[1]*Iudx1
    count_Iudx2_to_R = rho_r[2]*Iudx2
    count_Iudx3_to_R = rho_r[3]*Iudx3
    count_Iudx4_to_R = rho_r[4]*Iudx4
    count_Iudx5_to_R = rho_r[5]*Iudx5
    

    return(list(c(dS1, dS2, dS3, dS4, dS5, 
                  dE1, dE2, dE3, dE4, dE5,
                  dIudx1, dIudx2, dIudx3, dIudx4, dIudx4,
                  dIontime1, dIontime2, dIontime3, dIontime4, dIontime5,
                  dIlate1, dIlate2, dIlate3, dIlate4, dIlate5,
                  dInonad1, dInonad2, dInonad3, dInonad4, dInonad5,
                  dInoav1, dInoav2, dInoav3, dInoav4, dInoav5,
                  dR1, dR2, dR3, dR4, dR5,
                  count_Inoav1_to_R, count_Inoav2_to_R, count_Inoav3_to_R, count_Inoav4_to_R, count_Inoav5_to_R,
                  count_Ilate1_to_R, count_Ilate2_to_R, count_Ilate3_to_R, count_Ilate4_to_R, count_Ilate5_to_R,
                  count_Inonad1_to_R, count_Inonad2_to_R, count_Inonad3_to_R, count_Inonad4_to_R, count_Inonad5_to_R,
                  count_Iontime1_to_R, count_Iontime2_to_R, count_Iontime3_to_R, count_Iontime4_to_R, count_Iontime5_to_R,
                  count_Iudx1_to_R, count_Iudx2_to_R, count_Iudx3_to_R, count_Iudx4_to_R, count_Iudx5_to_R)))
  })
}
