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
# SET DURATIONS IN EACH COMPARTMENT - FIXED 
##############################################################################
#rate of exposed to infected (1/ latent period)
latent_pd = 1.5
gamma = 1 / latent_pd #days

#rate from undiagnosed infected to recovered
udx_to_r_pd = 3

#rates from undiagnosed infectious to care seeking infectious
udx_to_care_pd = 2.5 #days infectious before seeking care

#rates from care seeking infected to recovered
noav_rec_pd = .5
av_rec_pd = .5


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
            Iudx1 = 1000, Iudx2 = 1000, Iudx3 = 1000, Iudx4 = 1000, Iudx5 = 1000,
            Iontime1 = 0, Iontime2 = 0, Iontime3 = 0, Iontime4 = 0, Iontime5 = 0, 
            Ilate1 = 0, Ilate2 = 0, Ilate3 = 0, Ilate4 = 0, Ilate5 = 0, 
            Inonad1 = 0, Inonad2 = 0, Inonad3 = 0, Inonad4 = 0, Inonad5 = 0, 
            Inoav1 = 0, Inoav2 = 0, Inoav3 = 0, Inoav4 = 0, Inoav5 = 0,
            R1 = 0, R2 = 0, R3 = 0, R4 = 0, R5 = 0,
            count_Inoav1_to_R = 0, count_Inoav2_to_R = 0, count_Inoav3_to_R = 0, count_Inoav4_to_R = 0, count_Inoav5_to_R = 0,
            count_Ilate1_to_R = 0, count_Ilate2_to_R = 0, count_Ilate3_to_R = 0, count_Ilate4_to_R = 0, count_Ilate5_to_R = 0,
            count_Inonad1_to_R = 0, count_Inonad2_to_R = 0, count_Inonad3_to_R = 0, count_Inonad4_to_R = 0, count_Inonad5_to_R = 0,
            count_Iontime1_to_R = 0, count_Iontime2_to_R = 0, count_Iontime3_to_R = 0, count_Iontime4_to_R = 0, count_Iontime5_to_R = 0,
            count_Iudx1_to_R = 0, count_Iudx2_to_R = 0, count_Iudx3_to_R = 0, count_Iudx4_to_R = 0, count_Iudx5_to_R = 0)

