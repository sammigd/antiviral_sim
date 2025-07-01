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
            Iavzero1 = 0, Iavzero2 = 0, Iavzero3 = 0, Iavzero4 = 0 , Iavzero5 = 0,
            Iavone1 = 0, Iavone2 = 0, Iavone3 = 0, Iavone4 = 0, Iavone5 = 0,
            Iavtwo1 = 0, Iavtwo2 = 0, Iavtwo3 = 0, Iavtwo4 = 0, Iavtwo5 = 0,
            Iavthree1 = 0, Iavthree2 = 0, Iavthree3 = 0, Iavthree4 = 0, Iavthree5 = 0,
            R1 = 0, R2 = 0, R3 = 0, R4 = 0, R5 = 0,
            count_Iavzero_to_R1 = 0, count_Iavzero_to_R2 = 0, count_Iavzero_to_R3 = 0, count_Iavzero_to_R4 = 0, count_Iavzero_to_R5 = 0, 
            count_Iavone_to_R1 = 0, count_Iavone_to_R2 = 0, count_Iavone_to_R3 = 0, count_Iavone_to_R4 = 0, count_Iavone_to_R5 = 0, 
            count_Iavtwo_to_R1 = 0, count_Iavtwo_to_R2 = 0, count_Iavtwo_to_R3 = 0, count_Iavtwo_to_R4 = 0, count_Iavtwo_to_R5 = 0, 
            count_Iavthree_to_R1 = 0, count_Iavthree_to_R2 = 0, count_Iavthree_to_R3 = 0, count_Iavthree_to_R4 = 0, count_Iavthree_to_R5 = 0, 
            count_Isym_udx_to_R1 = 0, count_Isym_udx_to_R2 = 0, count_Isym_udx_to_R3 = 0, count_Isym_udx_to_R4 = 0, count_Isym_udx_to_R5 = 0, 
            count_Iasym_udx_to_R1 = 0, count_Iasym_udx_to_R2 = 0, count_Iasym_udx_to_R3 = 0, count_Iasym_udx_to_R4 = 0, count_Iasym_udx_to_R5 = 0
)

