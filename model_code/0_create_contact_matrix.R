#install.packages('fastmap')
#install.packages('devtools')

library(devtools)
#install_github('HHS/ASPR-flumodels')
library(flumodels)

#install_github('hrbrmstr/cdcfluview')
library(cdcfluview)

#install_github('HHS/ASPR-flumodelsutil')
library(flumodelsutil)


#loading age contact matrix from sinead's code
setwd('flu-vaccine-model')
hd_frac <- 0.75
source("0_simulation_functions.R")
source("1_call_simulation.R")
source("2_get_inputs.R")
setwd('..')

save(file = 'contact_matrix.RData', contacts)
save(file = 'age_fracs.RData', age_fracs)
save(file = 'popsize.RData', popsize)


rm(list = ls())

