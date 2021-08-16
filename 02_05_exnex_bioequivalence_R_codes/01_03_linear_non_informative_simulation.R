#--------------------------------------------------------------------------------------------------------------------------#
# R code for generating power under multiple scenario for the linear 2x2 crossover model non informative                   #
# Created by: Jonathan Jaeger - Date: 2020-12-17                                                                           #
# Last modified by: Jonathan Jaeger - Date: 2021-01-11                                                                     #
#--------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------- Source function -----------------------------------------------------#

   source(file = paste(directory,
                       '01_01_linear_non_informative_function.R',
                       sep = ''))


#----------------------------------------- Scenario 1 - Mean both group identical -----------------------------------------#

   system.time({

      N_subjid_per_sequence <- 10
      alpha <- c(4, 4)
      beta <- 0
      sigma_s <- 0.25
      sigma_epsilon <- 0.25

      name_object <- paste('power_linear_non_informative',
                           'N', N_subjid_per_sequence,
                           'alpha', paste(alpha, collapse = '_'),
                           'beta', beta,
                           'sigma_s', sigma_s,
                           'sigma_epsilon', sigma_epsilon,
                           sep = '_')

      power_int <- simulation_linear_non_informative(R = 500,
                                                     alpha = alpha, beta = beta,
                                                     sigma_s = sigma_s, sigma_epsilon = sigma_epsilon,
                                                     N_subjid_per_sequence = N_subjid_per_sequence,
                                                     seed = 20210111)
      
      assign(x = name_object,
      		 value = power_int)

   })

   save(list = name_object,
        file = paste(directory, name_object, '.RData', sep = ''))
