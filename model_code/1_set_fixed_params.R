##############################################################################
# LOAD IN FIXED POPULATION PARAMETERS
##############################################################################

# GET POPULATION AGE DISTRIBUTION ################
load('inputs/age_fracs.RData') #created in script 0
load('inputs/popsize.RData') #total pop denominator


# GET POPULATION CONTACT MATRIX ##################
load('inputs/contact_matrix.RData') #created in script 0

# COLLAPSE 5-12 AND 13-17 INTO A SINGLE AGE GROUP #
age_fracs_combined = c(age_fracs['0-4'], age_fracs['5-12'] + age_fracs['13-17'], age_fracs[c('18-49', '50-64', '65+')])
names(age_fracs_combined)[2] <- '5-17'

#combine by col
contacts_intermediate = contacts
contacts_intermediate[,2] = contacts_intermediate[,2] + contacts_intermediate[,3]
contacts_intermediate = contacts_intermediate[,-3]

#repeat by row
contacts_intermediate[2,] = contacts_intermediate[2,] + contacts_intermediate[3,]
contacts_combined = contacts_intermediate[-3,]

#get pop size by age group
age_pop_size = popsize*age_fracs_combined


##############################################################################
# SET INITIAL CONDITIONS OF MODEL #
##############################################################################
times = seq(from = 1, to = 365) 

N = popsize #pop size
N_byage = unname(popsize*age_fracs_combined)
state0 <- c(S1 = N_byage[1]-1000,
            S2 = N_byage[2]-1000,
            S3 = N_byage[3]-1000,
            S4 = N_byage[4]-1000,
            S5 = N_byage[5]-1000,
            E1 = 0,
            E2 = 0, 
            E3 = 0,
            E4 = 0, 
            E5 = 0,
            Isym_udx1 = 500, Isym_udx2 = 500, Isym_udx3 = 500, Isym_udx4 = 500, Isym_udx5 = 500,
            Iasym_udx1 = 500, Iasym_udx2 = 500, Iasym_udx3 = 500, Iasym_udx4 = 500, Iasym_udx5 = 500,
            Itime_adav1 = 0, Itime_adav2 = 0, Itime_adav3 = 0, Itime_adav4 = 0, Itime_adav5 = 0,
            Itime_part1 = 0, Itime_part2 = 0, Itime_part3 = 0, Itime_part4 = 0, Itime_part5 = 0,
            Itime_noav1 = 0, Itime_noav2 = 0, Itime_noav3 = 0, Itime_noav4 = 0, Itime_noav5 = 0,
            Ilate_av1 = 0, Ilate_av2 = 0, Ilate_av3 = 0, Ilate_av4 = 0, Ilate_av5 = 0,
            Ilate_noav1 = 0, Ilate_noav2 = 0, Ilate_noav3 = 0, Ilate_noav4 = 0, Ilate_noav5 = 0,
            Inocare_noav1 = 0, Inocare_noav2 = 0, Inocare_noav3 = 0, Inocare_noav4 = 0, Inocare_noav5 = 0,
            R1 = 0, R2 = 0, R3 = 0, R4 = 0, R5 = 0,
            count_Itime_adav1_to_R = 0, count_Itime_adav2_to_R = 0, count_Itime_adav3_to_R = 0, count_Itime_adav4_to_R = 0, count_Itime_adav5_to_R = 0,
            count_Itime_part1_to_R = 0, count_Itime_part2_to_R = 0, count_Itime_part3_to_R = 0, count_Itime_part4_to_R = 0, count_Itime_part5_to_R = 0,
            count_Itime_noav1_to_R = 0, count_Itime_noav2_to_R = 0, count_Itime_noav3_to_R = 0, count_Itime_noav4_to_R = 0, count_Itime_noav5_to_R = 0,
            count_Ilate_av1_to_R = 0, count_Ilate_av2_to_R = 0, count_Ilate_av3_to_R = 0, count_Ilate_av4_to_R = 0, count_Ilate_av5_to_R = 0,
            count_Ilate_noav1_to_R = 0, count_Ilate_noav2_to_R = 0, count_Ilate_noav3_to_R = 0, count_Ilate_noav4_to_R = 0, count_Ilate_noav5_to_R = 0,
            count_Inocare_noav1_to_R = 0, count_Inocare_noav2_to_R = 0, count_Inocare_noav3_to_R = 0, count_Inocare_noav4_to_R = 0, count_Inocare_noav5_to_R = 0,
            count_Iasym_udx1_to_R = 0, count_Iasym_udx2_to_R = 0, count_Iasym_udx3_to_R = 0, count_Iasym_udx4_to_R = 0, count_Iasym_udx5_to_R = 0
)

