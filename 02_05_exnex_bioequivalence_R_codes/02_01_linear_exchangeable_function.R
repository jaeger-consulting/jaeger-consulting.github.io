#--------------------------------------------------------------------------------------------------------------------------#
# R code for the function for the computation of power under exchangeable linear 2x2 crossover model                       #
# Created by: Jonathan Jaeger - Date: 2020-12-17                                                                           #
# Last modified by: Jonathan Jaeger - Date: 2020-12-18                                                                     #
#--------------------------------------------------------------------------------------------------------------------------#


#-------------------------------------------------------- Packages --------------------------------------------------------#

   require(parallel)
   require(rstan)


#------------------------------------------------ Function for simulations ------------------------------------------------#

   simulation_linear_exchangeable <- function(R,
                                              alpha, m_mean_alpha, t_mean_alpha,
                                              shape_sigma_alpha, scale_sigma_alpha,
                                              beta,
                                              sigma_s, sigma_epsilon,
                                              N_subjid_per_sequence,
                                              seed)
   {


   #-------------------------------- Lower and upper threshold for average bioequivalence ---------------------------------#

      lb_threshold <- -log(1.25)
      ub_threshold <- log(1.25)


   #---------------------------- Number of subjects, observation and unique subject identifier ----------------------------#

      N_subjid <- 2 * N_subjid_per_sequence
      N_observation <- 2 * N_subjid
      unique_subjid <- 1:N_subjid


   #------------------- Set seed for the data generations and the sampling of the posterior distribution ------------------#

      set.seed(seed)
      set.seed_seq <- sample.int(n = seed,
                                 size = R)


   #--------------------------------- Creating base dataset to be used for the simulation ---------------------------------#

      data_int_base <- expand.grid(subjid = unique_subjid,
                                   period = c(-1, 1))

      data_int_base $ sequence <- ifelse(data_int_base $ subjid <= N_subjid_per_sequence,
                                         1,
                                         2)

      data_int_base $ treatment <- ifelse(((data_int_base $ sequence %in% 1) &
                                           (data_int_base $ period %in% -1)) |
                                          ((data_int_base $ sequence %in% 2) &
                                           (data_int_base $ period %in% 1)),
                                          1,
                                          2)

      data_int_base $ treatment_fact <- as.factor(data_int_base $ treatment)

      data_int_base $ AUC_mean <- alpha[data_int_base $ treatment] +
                                  beta * data_int_base $ period


   #---------------------------------------------- Compiling the stan model -----------------------------------------------#

      model_stan <- rstan::stan_model(file = paste(directory,
                                                   '02_02_linear_exchangeable.stan',
                                                   sep = ''),
                                      verbose = FALSE)


   #-------------------------------------------- Preparing parallel computing ---------------------------------------------#

      no_cores <- detectCores(logical = TRUE)
      cl <- makeCluster(spec = no_cores - 1)
      clusterEvalQ(cl = cl,
                   expr = { require(rstan) })
      clusterExport(cl = cl,
                    varlist = c('data_int_base',
                                'unique_subjid', 'N_subjid', 'N_observation',
                                'alpha', 'm_mean_alpha', 't_mean_alpha',
                                'shape_sigma_alpha', 'scale_sigma_alpha',
                                'beta',
                                'sigma_epsilon', 'sigma_s',
                                'lb_threshold', 'ub_threshold',
                                'model_stan',
                                'set.seed_seq'),
                    envir = environment())


   #---------------------------------------- Running data generation and analysis -----------------------------------------#

      power_int <- parSapply(cl = cl,
                             X = 1:R,
                             FUN = { function(x)
                             {


                             #---------------------------------- Generating observations ----------------------------------#

                                set.seed(set.seed_seq[x])
                                data_int <- data_int_base
                                data_int <- merge(x = data_int,
                                                  y = data.frame(subjid = unique_subjid,
                                                                 s = rnorm(n = N_subjid,
                                                                           mean = 0.0,
                                                                           sd = sigma_s)),
                                                  by = 'subjid',
                                                  all.x = TRUE)
                                data_int $ AUC_obs <- data_int $ AUC_mean +
                                                      data_int $ s +
                                                      rnorm(n = N_observation,
                                                            mean = 0.0,
                                                            sd = sigma_epsilon)


                             #--------------------------- Data and inits for the stan sampling ----------------------------#

                                data_list <- list(N = N_observation,
                                                  N_formulation = 2,
                                                  N_subjid = N_subjid,
                                                  y = data_int $ AUC_obs,
                                                  index_formulation = data_int $ treatment,
                                                  period_covariate = data_int $ period,
                                                  index_subjid = data_int $ subjid,
                                                  m_mean_mu = m_mean_alpha,
                                                  t_mean_mu = t_mean_alpha,
                                                  shape_sigma_mu = shape_sigma_alpha,
                                                  scale_sigma_mu = scale_sigma_alpha)
                                inits_list <- function() { list(mu_tilde = 0 * alpha,
                                                                mean_mu = m_mean_alpha,
                                                                sigma_mu = 1 / (scale_sigma_alpha * (shape_sigma_alpha - 1)),
                                                                beta = beta,
                                                                sigma_epsilon = sigma_epsilon,
                                                                sigma_s = sigma_s,
                                                                s = rep(x = 0, times = N_subjid)) }


                             #--------------------------- Sample of the posterior distribution ----------------------------#

                                fit <- sampling(object = model_stan,
                                                data = data_list,
                                                init = inits_list,
                                                pars = c('mu', 'mean_mu', 'sigma_mu', 'beta', 'sigma_epsilon', 'sigma_s'),
                                                seed = set.seed_seq[x],
                                                chains = 1,
                                                iter = 22000,
                                                warmup = 2000)
                                fit_sample <- as.matrix(fit)


                             #---------------------------- Checking for average bioequivalence ----------------------------#

                                delta <- apply(X = fit_sample[, c('mu[1]', 'mu[2]')], MARGIN = 1, diff)
                                CI_delta <- quantile(x = delta, probs = c(0.05, 0.95))
                                average_bioequivalence <- (CI_delta[1] >= lb_threshold) & (CI_delta[2] <= ub_threshold)


                             #---------------------- Checking coverage mu[1] and mu[2] versus alpha -----------------------#

                                CI_mu <- apply(X = fit_sample[, c('mu[1]', 'mu[2]')],
                                               MARGIN = 2,
                                               FUN = { function(x) quantile(x = x, probs = c(0.05, 0.95)) })
                                coverage_mu_1 <- (CI_mu[1, 1] <= alpha[1]) & (CI_mu[2, 1] >= alpha[1])
                                coverage_mu_2 <- (CI_mu[1, 2] <= alpha[2]) & (CI_mu[2, 2] >= alpha[2])


                             #--------------------------------- Return relevant elements ----------------------------------#

                                return(c(average_bioequivalence, coverage_mu_1, coverage_mu_2))
                             } })


   #--------------------------------------------- Stopping parallel computing ---------------------------------------------#

      stopCluster(cl = cl)


   #--------------------------------------- Return the vector of the study results ----------------------------------------#

      return(power_int)
   }
