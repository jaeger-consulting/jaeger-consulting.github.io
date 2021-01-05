#--------------------------------------------------------------------------------------------------------------------------#
# R code for generating power under multiple scenario for the linear 2x2 crossover model exchangeable                      #
# Created by: Jonathan Jaeger - Date: 2020-12-17                                                                           #
# Last modified by: Jonathan Jaeger - Date: 2020-12-18                                                                     #
#--------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------- Source function -----------------------------------------------------#

   source(file = paste(directory,
                       '02_01_linear_exchangeable_function.R',
                       sep = ''))


#----------------------------------------- Scenario 1 - Mean both group identical -----------------------------------------#

   system.time({

      N_subjid_per_sequence <- 12
      alpha <- c(4, 6)
      m_mean_alpha <- 4
      t_mean_alpha <- 0.1
      shape_sigma_alpha <- 10 #11.5 #18
      rate_sigma_alpha <- 2 #2 #3.4
      scale_sigma_alpha <- 1 / rate_sigma_alpha
      beta <- 0
      sigma_s <- 0.25
      sigma_epsilon <- 0.25

      name_object <- paste('power_linear_exchangeable',
                           'N', N_subjid_per_sequence,
                           'alpha', paste(alpha, collapse = '_'),
                           'beta', beta,
                           'sigma_s', sigma_s,
                           'sigma_epsilon', sigma_epsilon,
                           sep = '_')

      power_int <- simulation_linear_exchangeable(R = 50,
                                                  alpha = alpha,
                                                  m_mean_alpha = m_mean_alpha, t_mean_alpha = t_mean_alpha,
                                                  shape_sigma_alpha = shape_sigma_alpha, scale_sigma_alpha = scale_sigma_alpha,
                                                  beta = beta,
                                                  sigma_s = sigma_s, sigma_epsilon = sigma_epsilon,
                                                  N_subjid_per_sequence = N_subjid_per_sequence,
                                                  seed = 20201217)

      assign(x = name_object,
      		 value = power_int)

   })

   save(list = name_object,
        file = paste(directory, name_object, '.RData', sep = ''))
