##################################################
# WRITE MODEL FUNCTION ###########################
##################################################

flu.model <- function(t, state, parameters){
  with(as.list(c(state,parameters)),{
    
    #add seasonality to betas
    betas_noav = betas_noav * (1 + .1*sin(2*pi*(t/365)))
    betas_av   = betas_av   * (1 + .1*sin(2*pi*(t/365)))
    betas_asym = 0.57*betas_noav
    betas_three = betas_noav

    #TOTAL POPULATION - SPLIT INTO 5 AGE GROUPS
    N <- S1 + E1 + Isym_udx1 + Iasym_udx1 + Iavzero1 + Iavone1 + Iavtwo1 + Iavthree1 + R1 +
         S2 + E2 + Isym_udx2 + Iasym_udx2 + Iavzero2 + Iavone2 + Iavtwo2 + Iavthree2 + R2 +
         S3 + E3 + Isym_udx3 + Iasym_udx3 + Iavzero3 + Iavone3 + Iavtwo3 + Iavthree3 + R3 +
         S4 + E4 + Isym_udx4 + Iasym_udx4 + Iavzero4 + Iavone4 + Iavtwo4 + Iavthree4 + R4 +
         S5 + E5 + Isym_udx5 + Iasym_udx5 + Iavzero5 + Iavone5 + Iavtwo5 + Iavthree5 + R5
    
    vec_Isym_noav = c(Isym_udx1, Isym_udx2, Isym_udx3, Isym_udx4, Isym_udx5)
    
    vec_Isym_av = c(Iavzero1 + Iavone1 + Iavtwo1 ,
                    Iavzero2 + Iavone2 + Iavtwo2 ,
                    Iavzero3 + Iavone3 + Iavtwo3 ,
                    Iavzero4 + Iavone4 + Iavtwo4 ,
                    Iavzero5 + Iavone5 + Iavtwo5 )
    
    vec_Iasym = c(Iasym_udx1, Iasym_udx2, Iasym_udx3, Iasym_udx4, Iasym_udx5)
    
    vec_Isym_three = c(Iavthree1, Iavthree2, Iavthree3, Iavthree4, Iavthree5) #special case: rcvd av on day three after symptom onset

    dS1 <- -sum(betas_noav[1,]*S1*vec_Isym_noav) - sum(betas_asym[1,]*S1*vec_Iasym) - sum(betas_av[1,]*S1*vec_Isym_av) - sum(betas_three[1,]*S1*vec_Isym_three) 
    dS2 <- -sum(betas_noav[2,]*S2*vec_Isym_noav) - sum(betas_asym[2,]*S2*vec_Iasym) - sum(betas_av[2,]*S2*vec_Isym_av) - sum(betas_three[2,]*S2*vec_Isym_three)
    dS3 <- -sum(betas_noav[3,]*S3*vec_Isym_noav) - sum(betas_asym[3,]*S3*vec_Iasym) - sum(betas_av[3,]*S3*vec_Isym_av) - sum(betas_three[3,]*S3*vec_Isym_three) 
    dS4 <- -sum(betas_noav[4,]*S4*vec_Isym_noav) - sum(betas_asym[4,]*S4*vec_Iasym) - sum(betas_av[4,]*S4*vec_Isym_av) - sum(betas_three[4,]*S4*vec_Isym_three)
    dS5 <- -sum(betas_noav[5,]*S5*vec_Isym_noav) - sum(betas_asym[5,]*S5*vec_Iasym) - sum(betas_av[5,]*S5*vec_Isym_av) - sum(betas_three[5,]*S5*vec_Isym_three)
    
    dE1 <-  sum(betas_noav[1,]*S1*vec_Isym_noav) + sum(betas_asym[1,]*S1*vec_Iasym) + sum(betas_av[1,]*S1*vec_Isym_av) + sum(betas_three[1,]*S1*vec_Isym_three) - gamma_asym*E1 - gamma_sym*E1
    dE2 <-  sum(betas_noav[2,]*S2*vec_Isym_noav) + sum(betas_asym[2,]*S2*vec_Iasym) + sum(betas_av[2,]*S2*vec_Isym_av) + sum(betas_three[2,]*S2*vec_Isym_three) - gamma_asym*E2 - gamma_sym*E2
    dE3 <-  sum(betas_noav[3,]*S3*vec_Isym_noav) + sum(betas_asym[3,]*S3*vec_Iasym) + sum(betas_av[3,]*S3*vec_Isym_av) + sum(betas_three[3,]*S3*vec_Isym_three) - gamma_asym*E3 - gamma_sym*E3
    dE4 <-  sum(betas_noav[4,]*S4*vec_Isym_noav) + sum(betas_asym[4,]*S4*vec_Iasym) + sum(betas_av[4,]*S4*vec_Isym_av) + sum(betas_three[4,]*S4*vec_Isym_three) - gamma_asym*E4 - gamma_sym*E4
    dE5 <-  sum(betas_noav[5,]*S5*vec_Isym_noav) + sum(betas_asym[5,]*S5*vec_Iasym) + sum(betas_av[5,]*S5*vec_Isym_av) + sum(betas_three[5,]*S5*vec_Isym_three) - gamma_asym*E5 - gamma_sym*E5
    
    dIasym_udx1 = gamma_asym*E1 - rho_asym_r*Iasym_udx1                         # undiagnosed asymptomatic compartments
    dIasym_udx2 = gamma_asym*E2 - rho_asym_r*Iasym_udx2 
    dIasym_udx3 = gamma_asym*E3 - rho_asym_r*Iasym_udx3
    dIasym_udx4 = gamma_asym*E4 - rho_asym_r*Iasym_udx4 
    dIasym_udx5 = gamma_asym*E5 - rho_asym_r*Iasym_udx5 
    
    dIsym_udx1 <- gamma_sym*E1 - Isym_udx1*(rho_avzero[1] + rho_avone[1] + rho_avtwo[1] + rho_avthree[1] + rho_r[1])      # udx symptomatic compartments
    dIsym_udx2 <- gamma_sym*E2 - Isym_udx2*(rho_avzero[2] + rho_avone[2] + rho_avtwo[2] + rho_avthree[2] + rho_r[2])      # from: E
    dIsym_udx3 <- gamma_sym*E3 - Isym_udx3*(rho_avzero[3] + rho_avone[3] + rho_avtwo[3] + rho_avthree[3] + rho_r[3])      # to: av on days 0-3, R
    dIsym_udx4 <- gamma_sym*E4 - Isym_udx4*(rho_avzero[4] + rho_avone[4] + rho_avtwo[4] + rho_avthree[4] + rho_r[4])      #     
    dIsym_udx5 <- gamma_sym*E5 - Isym_udx5*(rho_avzero[5] + rho_avone[5] + rho_avtwo[5] + rho_avthree[5] + rho_r[5])      #     
    
    dIavzero1 <- rho_avzero[1]*Isym_udx1 - mu_avzero[1]*Iavzero1
    dIavzero2 <- rho_avzero[2]*Isym_udx2 - mu_avzero[2]*Iavzero2
    dIavzero3 <- rho_avzero[3]*Isym_udx3 - mu_avzero[3]*Iavzero3
    dIavzero4 <- rho_avzero[4]*Isym_udx4 - mu_avzero[4]*Iavzero4
    dIavzero5 <- rho_avzero[5]*Isym_udx5 - mu_avzero[5]*Iavzero5
    
    dIavone1 <- rho_avone[1]*Isym_udx1 - mu_avone[1]*Iavone1
    dIavone2 <- rho_avone[2]*Isym_udx2 - mu_avone[2]*Iavone2
    dIavone3 <- rho_avone[3]*Isym_udx3 - mu_avone[3]*Iavone3
    dIavone4 <- rho_avone[4]*Isym_udx4 - mu_avone[4]*Iavone4
    dIavone5 <- rho_avone[5]*Isym_udx5 - mu_avone[5]*Iavone5
    
    dIavtwo1 <- rho_avtwo[1]*Isym_udx1 - mu_avtwo[1]*Iavtwo1
    dIavtwo2 <- rho_avtwo[2]*Isym_udx2 - mu_avtwo[2]*Iavtwo2
    dIavtwo3 <- rho_avtwo[3]*Isym_udx3 - mu_avtwo[3]*Iavtwo3
    dIavtwo4 <- rho_avtwo[4]*Isym_udx4 - mu_avtwo[4]*Iavtwo4
    dIavtwo5 <- rho_avtwo[5]*Isym_udx5 - mu_avtwo[5]*Iavtwo5
    
    dIavthree1 <- rho_avthree[1]*Isym_udx1 - mu_avthree[1]*Iavthree1
    dIavthree2 <- rho_avthree[2]*Isym_udx2 - mu_avthree[2]*Iavthree2
    dIavthree3 <- rho_avthree[3]*Isym_udx3 - mu_avthree[3]*Iavthree3
    dIavthree4 <- rho_avthree[4]*Isym_udx4 - mu_avthree[4]*Iavthree4
    dIavthree5 <- rho_avthree[5]*Isym_udx5 - mu_avthree[5]*Iavthree5
    
    dR1 <- mu_avzero[1]*Iavzero1 + mu_avone[1]*Iavone1 + mu_avtwo[1]*Iavtwo1 + mu_avthree[1]*Iavthree1 + rho_r[1]*Isym_udx1 + rho_asym_r*Iasym_udx1
    dR2 <- mu_avzero[2]*Iavzero2 + mu_avone[2]*Iavone2 + mu_avtwo[2]*Iavtwo2 + mu_avthree[2]*Iavthree2 + rho_r[2]*Isym_udx2 + rho_asym_r*Iasym_udx2
    dR3 <- mu_avzero[3]*Iavzero3 + mu_avone[3]*Iavone3 + mu_avtwo[3]*Iavtwo3 + mu_avthree[3]*Iavthree3 + rho_r[3]*Isym_udx3 + rho_asym_r*Iasym_udx3
    dR4 <- mu_avzero[4]*Iavzero4 + mu_avone[4]*Iavone4 + mu_avtwo[4]*Iavtwo4 + mu_avthree[4]*Iavthree4 + rho_r[4]*Isym_udx4 + rho_asym_r*Iasym_udx4
    dR5 <- mu_avzero[5]*Iavzero5 + mu_avone[5]*Iavone5 + mu_avtwo[5]*Iavtwo5 + mu_avthree[5]*Iavthree5 + rho_r[5]*Isym_udx5 + rho_asym_r*Iasym_udx5
    
  
    #COUNT TRANSITIONS FROM EACH COMPARTMENT TO RECOVERED
    count_Iavzero_to_R1 = mu_avzero[1]*Iavzero1 
    count_Iavzero_to_R2 = mu_avzero[2]*Iavzero2
    count_Iavzero_to_R3 = mu_avzero[3]*Iavzero3
    count_Iavzero_to_R4 = mu_avzero[4]*Iavzero4
    count_Iavzero_to_R5 = mu_avzero[5]*Iavzero5
    
    count_Iavone_to_R1 = mu_avone[1]*Iavone1 
    count_Iavone_to_R2 = mu_avone[2]*Iavone2
    count_Iavone_to_R3 = mu_avone[3]*Iavone3
    count_Iavone_to_R4 = mu_avone[4]*Iavone4
    count_Iavone_to_R5 = mu_avone[5]*Iavone5
    
    count_Iavtwo_to_R1 = mu_avtwo[1]*Iavtwo1 
    count_Iavtwo_to_R2 = mu_avtwo[2]*Iavtwo2
    count_Iavtwo_to_R3 = mu_avtwo[3]*Iavtwo3
    count_Iavtwo_to_R4 = mu_avtwo[4]*Iavtwo4
    count_Iavtwo_to_R5 = mu_avtwo[5]*Iavtwo5
    
    count_Iavthree_to_R1 = mu_avthree[1]*Iavthree1 
    count_Iavthree_to_R2 = mu_avthree[2]*Iavthree2
    count_Iavthree_to_R3 = mu_avthree[3]*Iavthree3
    count_Iavthree_to_R4 = mu_avthree[4]*Iavthree4
    count_Iavthree_to_R5 = mu_avthree[5]*Iavthree5
    
    count_Isym_udx_to_R1 = rho_r[1]*Isym_udx1
    count_Isym_udx_to_R2 = rho_r[2]*Isym_udx2
    count_Isym_udx_to_R3 = rho_r[3]*Isym_udx3
    count_Isym_udx_to_R4 = rho_r[4]*Isym_udx4
    count_Isym_udx_to_R5 = rho_r[5]*Isym_udx5
    
    count_Iasym_udx_to_R1 = rho_asym_r*Iasym_udx1
    count_Iasym_udx_to_R2 = rho_asym_r*Iasym_udx2
    count_Iasym_udx_to_R3 = rho_asym_r*Iasym_udx3
    count_Iasym_udx_to_R4 = rho_asym_r*Iasym_udx4
    count_Iasym_udx_to_R5 = rho_asym_r*Iasym_udx5

    return(list(c(dS1, dS2, dS3, dS4, dS5, 
                  dE1, dE2, dE3, dE4, dE5,
                  dIsym_udx1, dIsym_udx2, dIsym_udx3, dIsym_udx4, dIsym_udx5,
                  dIasym_udx1, dIasym_udx2, dIasym_udx3, dIasym_udx4, dIasym_udx5,
                  dIavzero1, dIavzero2, dIavzero3, dIavzero4, dIavzero5,
                  dIavone1, dIavone2, dIavone3, dIavone4, dIavone5,
                  dIavtwo1, dIavtwo2, dIavtwo3, dIavtwo4, dIavtwo5,
                  dIavthree1, dIavthree2, dIavthree3, dIavthree4, dIavthree5,
                  dR1, dR2, dR3, dR4, dR5,
                  count_Iavzero_to_R1, count_Iavzero_to_R2, count_Iavzero_to_R3, count_Iavzero_to_R4, count_Iavzero_to_R5,
                  count_Iavone_to_R1, count_Iavone_to_R2, count_Iavone_to_R3, count_Iavone_to_R4, count_Iavone_to_R5,
                  count_Iavtwo_to_R1, count_Iavtwo_to_R2, count_Iavtwo_to_R3, count_Iavtwo_to_R4, count_Iavtwo_to_R5,
                  count_Iavthree_to_R1, count_Iavthree_to_R2, count_Iavthree_to_R3, count_Iavthree_to_R4, count_Iavthree_to_R5,
                  count_Isym_udx_to_R1, count_Isym_udx_to_R2, count_Isym_udx_to_R3, count_Isym_udx_to_R4, count_Isym_udx_to_R5,
                  count_Iasym_udx_to_R1, count_Iasym_udx_to_R2, count_Iasym_udx_to_R3, count_Iasym_udx_to_R4, count_Iasym_udx_to_R5
                  )))
  })
}
