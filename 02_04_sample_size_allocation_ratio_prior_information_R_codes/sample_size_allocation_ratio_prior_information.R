#--------------------------------------------------------------------------------------------------------------------------#
# R code for sample size calculation in Bayesian setting in presence of prior information                                  #
# Gaussian data distribution                                                                                               #
# Created by: Jonathan Jaeger - Date: 2020-10-27                                                                           #
# Last modified by: Jonathan Jaeger - Date: 2020-10-28                                                                     #
#--------------------------------------------------------------------------------------------------------------------------#


#-------------------------------------------------------- Packages --------------------------------------------------------#

   require(parallel)
   require(StanHeaders)
   require(rstan)


#------------------------------------------------ Compiling the stan model ------------------------------------------------#

   model_stan <- rstan::stan_model(file = paste('G:/03_Projects/001_Website/jaeger-consulting.github.io/documents/',
                                                'sample_size_allocation_ratio_prior_information.stan',
                                                sep = ''),
                                   verbose = FALSE)


#----------------------------------------- Common assumptions accross simulations -----------------------------------------#


   #---------------------------------------------------- Control group ----------------------------------------------------#

      mu_ctl <- 5
      sigma_ctl <- 2

      mean_mu_ctl_non_informative <- 0
      mean_mu_ctl_informative <- mu_ctl
      sd_mu_ctl_non_informative <- 1e3
      sd_mu_ctl_informative <- 1.0


   #------------------------------------------------- Experimental group --------------------------------------------------#

      mu_exp <- 7
      sigma_exp <- 2

      mean_mu_exp_non_informative <- 0
      sd_mu_exp_non_informative <- 1e3


   #----------------------------------------- Target operational characteristics ------------------------------------------#

      superiority_threshold <- 0.975
      power_target <- 0.8


   #------------------------------------------------ Number of replicates -------------------------------------------------#

      R <- 2e4


   #-------------------------------------- Number of iterations and warmup for MCMC ---------------------------------------#

      iter <- 6e3
      warmup <- 1e3


#--------------------------------------- Sample size for target power in 1:1 ratio ----------------------------------------#


   #-------------------------------------------- Setting up parallel computing --------------------------------------------#

      no_cores <- detectCores(logical = TRUE)
      cl <- makeCluster(no_cores - 2)
      clusterExport(cl = cl,
                    varlist = c('model_stan',
                                'mean_mu_ctl_non_informative', 'mean_mu_exp_non_informative',
                                'sd_mu_ctl_non_informative', 'sd_mu_exp_non_informative',
                                'superiority_threshold', 'power_target',
                                'iter', 'warmup'))
      clusterEvalQ(cl = cl,
                   expr = {require(StanHeaders); require(rstan)})


   #--------------------------------------------- Sample size starting values ---------------------------------------------#

      n_ctl <- 20
      n_exp <- n_ctl


   #------------------------------------------------- Storage of results --------------------------------------------------#

      mean_sim_result_1 <- NULL


   #---------------------------------------- Perform simulation until target power ----------------------------------------#

      sim_result <- 1

      while(mean(sim_result > superiority_threshold) >= power_target)
      {
         n_ctl <- n_ctl - 1
         n_exp <- n_exp - 1

         cat('Scenario 1 - n_ctl and n_exp:', n_ctl, '\n')
         flush.console()

         scenario_mat <- cbind(n_ctl = rep(x = n_ctl, R),
                               mu_ctl = rep(x = mu_ctl, R),
                               sigma_ctl = rep(x = sigma_ctl, R),
                               n_exp = rep(x = n_exp, R),
                               mu_exp = rep(x = mu_exp, R),
                               sigma_exp = rep(x = sigma_exp, R))

         sim_result <- parApply(cl = cl,
                                X = scenario_mat,
                                MARGIN = 1,
                                FUN = { function(x)
                                        {
                                           # control
                                           y_ctl <- rnorm(n = x[1], mean = x[2], sd = x[3])
                                           fit_ctl <- sampling(object = model_stan,
                                                               data = list(N = x[1],
                                                                           y = as.numeric(y_ctl),
                                                                           mean_mu = mean_mu_ctl_non_informative,
                                                                           sd_mu = sd_mu_ctl_non_informative),
                                                               init = function() { list(mu = mean(y_ctl),
                                                                                        tau = 1 / var(y_ctl)) },
                                                               pars = c('mu', 'sigma'),
                                                               iter = iter,
                                                               warmup = warmup,
                                                               chains = 1,
                                                               cores = 1)
                                           fit_ctl <- as.matrix(fit_ctl)[, 'mu']

                                           # experimental
                                           y_exp <- rnorm(n = x[4], mean = x[5], sd = x[6])
                                           fit_exp <- sampling(object = model_stan,
                                                               data = list(N = x[4],
                                                                           y = as.numeric(y_exp),
                                                                           mean_mu = mean_mu_exp_non_informative,
                                                                           sd_mu = sd_mu_exp_non_informative),
                                                               init = function() { list(mu = mean(y_exp),
                                                                                        tau = 1 / var(y_exp)) },
                                                               pars = c('mu', 'sigma'),
                                                               iter = iter,
                                                               warmup = warmup,
                                                               chains = 1,
                                                               cores = 1)
                                           fit_exp <- as.matrix(fit_exp)[, 'mu']

                                           # comparison control vs experimental
                                           p_ctl_exp <- mean(fit_ctl < fit_exp)

                                        }})
         mean_sim_result_1 <- cbind(mean_sim_result_1,
                                    c(n_ctl, n_exp, mean(sim_result > superiority_threshold)))

         if (R >= 1e4) { Sys.sleep(time = 60) }
      }


   #-------------------------------------------- Closing up parallel computing --------------------------------------------#

      stopCluster(cl = cl)


#--------- Power with varying allocation ratio and no prior information for fixed sample size computed previously ---------#


   #-------------------------------------------- Setting up parallel computing --------------------------------------------#

      no_cores <- detectCores(logical = TRUE)
      cl <- makeCluster(no_cores - 2)
      clusterExport(cl = cl,
                    varlist = c('model_stan',
                                'mean_mu_ctl_non_informative', 'mean_mu_exp_non_informative',
                                'sd_mu_ctl_non_informative', 'sd_mu_exp_non_informative',
                                'superiority_threshold', 'power_target',
                                'iter', 'warmup'))
      clusterEvalQ(cl = cl,
                   expr = {require(StanHeaders); require(rstan)})


   #------------------------------------------------- Overall sample size -------------------------------------------------#

      n_ctl <- mean_sim_result_1[1, max(which(mean_sim_result_1[3, ] >= power_target))]
      n_exp <- mean_sim_result_1[2, max(which(mean_sim_result_1[3, ] >= power_target))]
      n_total <- n_ctl + n_exp


   #------------------------------------------------- Storage of results --------------------------------------------------#

      mean_sim_result_2 <- NULL


   #---------------------------------- Perform simulation over range of allocation ratio ----------------------------------#

      range_n_ctl <- (ceiling(0.5 * n_total) + (-floor(0.5 * (0.5 * n_total - 1)):floor(0.5 * (0.5 * n_total - 1))))

      for (n_ctl in range_n_ctl)
      {
         cat('Scenario 2 - n_ctl: ', n_ctl, '\n')
         flush.console()

         n_exp <- n_total - n_ctl

         scenario_mat <- cbind(n_ctl = rep(x = n_ctl, R),
                               mu_ctl = rep(x = mu_ctl, R),
                               sigma_ctl = rep(x = sigma_ctl, R),
                               n_exp = rep(x = n_exp, R),
                               mu_exp = rep(x = mu_exp, R),
                               sigma_exp = rep(x = sigma_exp, R))

         sim_result <- parApply(cl = cl,
                                X = scenario_mat,
                                MARGIN = 1,
                                FUN = { function(x)
                                        {
                                           # control
                                           y_ctl <- rnorm(n = x[1], mean = x[2], sd = x[3])
                                           fit_ctl <- sampling(object = model_stan,
                                                               data = list(N = x[1],
                                                                           y = as.numeric(y_ctl),
                                                                           mean_mu = mean_mu_ctl_non_informative,
                                                                           sd_mu = sd_mu_ctl_non_informative),
                                                               init = function() { list(mu = mean(y_ctl),
                                                                                        tau = 1 / var(y_ctl)) },
                                                               pars = c('mu', 'sigma'),
                                                               iter = iter,
                                                               warmup = warmup,
                                                               chains = 1,
                                                               cores = 1)
                                           fit_ctl <- as.matrix(fit_ctl)[, 'mu']

                                           # experimental
                                           y_exp <- rnorm(n = x[4], mean = x[5], sd = x[6])
                                           fit_exp <- sampling(object = model_stan,
                                                               data = list(N = x[4],
                                                                           y = as.numeric(y_exp),
                                                                           mean_mu = mean_mu_exp_non_informative,
                                                                           sd_mu = sd_mu_exp_non_informative),
                                                               init = function() { list(mu = mean(y_exp),
                                                                                        tau = 1 / var(y_exp)) },
                                                               pars = c('mu', 'sigma'),
                                                               iter = iter,
                                                               warmup = warmup,
                                                               chains = 1,
                                                               cores = 1)
                                           fit_exp <- as.matrix(fit_exp)[, 'mu']

                                           # comparison control vs experimental
                                           p_ctl_exp <- mean(fit_ctl < fit_exp)

                                        }})

         mean_sim_result_2 <- cbind(mean_sim_result_2,
                                    c(n_ctl, n_exp, mean(sim_result > superiority_threshold)))

         if (R >= 1e4) { Sys.sleep(time = 60) }
      }


   #-------------------------------------------- Closing up parallel computing --------------------------------------------#

      stopCluster(cl = cl)


#-------- Power with varying allocation ratio and with prior information for fixed sample size computed previously --------#


   #-------------------------------------------- Setting up parallel computing --------------------------------------------#

      no_cores <- detectCores(logical = TRUE)
      cl <- makeCluster(no_cores - 2)
      clusterExport(cl = cl,
                    varlist = c('model_stan',
                                'mean_mu_ctl_informative', 'mean_mu_exp_non_informative',
                                'sd_mu_ctl_informative', 'sd_mu_exp_non_informative',
                                'superiority_threshold', 'power_target',
                                'iter', 'warmup'))
      clusterEvalQ(cl = cl,
                   expr = {require(StanHeaders); require(rstan)})


   #------------------------------------------------- Storage of results --------------------------------------------------#

      mean_sim_result_3 <- NULL


   #---------------------------------- Perform simulation over range of allocation ratio ----------------------------------#

      for (n_ctl in range_n_ctl)
      {
         cat('Scenario 3 - n_ctl: ', n_ctl, '\n')
         flush.console()

         n_exp <- n_total - n_ctl

         scenario_mat <- cbind(n_ctl = rep(x = n_ctl, R),
                               mu_ctl = rep(x = mu_ctl, R),
                               sigma_ctl = rep(x = sigma_ctl, R),
                               n_exp = rep(x = n_exp, R),
                               mu_exp = rep(x = mu_exp, R),
                               sigma_exp = rep(x = sigma_exp, R))

         sim_result <- parApply(cl = cl,
                                X = scenario_mat,
                                MARGIN = 1,
                                FUN = { function(x)
                                        {
                                           # control
                                           y_ctl <- rnorm(n = x[1], mean = x[2], sd = x[3])
                                           fit_ctl <- sampling(object = model_stan,
                                                               data = list(N = x[1],
                                                                           y = as.numeric(y_ctl),
                                                                           mean_mu = mean_mu_ctl_informative,
                                                                           sd_mu = sd_mu_ctl_informative),
                                                               init = function() { list(mu = mean(y_ctl),
                                                                                        tau = 1 / var(y_ctl)) } ,
                                                               pars = c('mu', 'sigma'),
                                                               iter = iter,
                                                               warmup = warmup,
                                                               chains = 1,
                                                               cores = 1)
                                           fit_ctl <- as.matrix(fit_ctl)[, 'mu']

                                           # experimental
                                           y_exp <- rnorm(n = x[4], mean = x[5], sd = x[6])
                                           fit_exp <- sampling(object = model_stan,
                                                               data = list(N = x[4],
                                                                           y = as.numeric(y_exp),
                                                                           mean_mu = mean_mu_exp_non_informative,
                                                                           sd_mu = sd_mu_exp_non_informative),
                                                               init = function() { list(mu = mean(y_exp),
                                                                                        tau = 1 / var(y_exp)) } ,
                                                               pars = c('mu', 'sigma'),
                                                               iter = iter,
                                                               warmup = warmup,
                                                               chains = 1,
                                                               cores = 1)
                                           fit_exp <- as.matrix(fit_exp)[, 'mu']

                                           # comparison control vs experimental
                                           p_ctl_exp <- mean(fit_ctl < fit_exp)

                                        }})

         mean_sim_result_3 <- cbind(mean_sim_result_3,
                                      c(n_ctl, n_exp, mean(sim_result > superiority_threshold)))

         if (R >= 1e4) { Sys.sleep(time = 60) }

      }


   #-------------------------------------------- Closing up parallel computing --------------------------------------------#

      stopCluster(cl = cl)


#------------------- Sample size optimization with varying allocation ratio to target expected power ----------------------#


   #-------------------------------------------- Setting up parallel computing --------------------------------------------#

      no_cores <- detectCores(logical = TRUE)
      cl <- makeCluster(no_cores - 2)
      clusterExport(cl = cl,
                    varlist = c('model_stan',
                                'mean_mu_ctl_informative', 'mean_mu_exp_non_informative',
                                'sd_mu_ctl_informative', 'sd_mu_exp_non_informative',
                                'superiority_threshold', 'power_target',
                                'iter', 'warmup'))
      clusterEvalQ(cl = cl,
                   expr = {require(StanHeaders); require(rstan)})


   #----------------------------------------- Overall sample size starting values -----------------------------------------#

      n_total <- 30


   #------------------------------------------------- Storage of results --------------------------------------------------#

      j <- 1
      mean_sim_result_4 <- list()

      mean_sim_result_4[[j]] <- matrix(c(NA, NA, 1), ncol = 1)


   #---------------------------------------- Perform simulation until target power ----------------------------------------#

      while (any(mean_sim_result_4[[j]][3, ] >= power_target))
      {
         j <- j + 1

         length(mean_sim_result_4) <- j

         n_total <- n_total - 1
         range_n_ctl <- (ceiling(0.5 * n_total) + (-floor(0.7 * (0.5 * n_total - 1)):0))

         for (n_ctl in range_n_ctl)
         {
            cat('Scenario 4 - n_total:', n_total, '- n_ctl:', n_ctl, '\n')
            flush.console()

            n_exp <- n_total - n_ctl

            scenario_mat <- cbind(n_ctl = rep(x = n_ctl, R),
                                  mu_ctl = rep(x = mu_ctl, R),
                                  sigma_ctl = rep(x = sigma_ctl, R),
                                  n_exp = rep(x = n_exp, R),
                                  mu_exp = rep(x = mu_exp, R),
                                  sigma_exp = rep(x = sigma_exp, R))

            sim_result <- parApply(cl = cl,
                                   X = scenario_mat,
                                   MARGIN = 1,
                                   FUN = { function(x)
                                           {
                                              # control
                                              y_ctl <- rnorm(n = x[1], mean = x[2], sd = x[3])
                                              fit_ctl <- sampling(object = model_stan,
                                                                  data = list(N = x[1],
                                                                              y = as.numeric(y_ctl),
                                                                              mean_mu = mean_mu_ctl_informative,
                                                                              sd_mu = sd_mu_ctl_informative),
                                                                  init = function() { list(mu = mean(y_ctl),
                                                                                           tau = 1 / var(y_ctl)) } ,
                                                                  pars = c('mu', 'sigma'),
                                                                  iter = iter,
                                                                  warmup = warmup,
                                                                  chains = 1,
                                                                  cores = 1)
                                              fit_ctl <- as.matrix(fit_ctl)[, 'mu']

                                              # experimental
                                              y_exp <- rnorm(n = x[4], mean = x[5], sd = x[6])
                                              fit_exp <- sampling(object = model_stan,
                                                                  data = list(N = x[4],
                                                                              y = as.numeric(y_exp),
                                                                              mean_mu = mean_mu_exp_non_informative,
                                                                              sd_mu = sd_mu_exp_non_informative),
                                                                  init = function() { list(mu = mean(y_exp),
                                                                                           tau = 1 / var(y_exp)) } ,
                                                                  pars = c('mu', 'sigma'),
                                                                  iter = iter,
                                                                  warmup = warmup,
                                                                  chains = 1,
                                                                  cores = 1)
                                              fit_exp <- as.matrix(fit_exp)[, 'mu']

                                              # comparison control vs experimental
                                              p_ctl_exp <- mean(fit_ctl < fit_exp)

                                           }})

            mean_sim_result_4[[j]] <- cbind(mean_sim_result_4[[j]],
                                            c(n_ctl, n_exp, mean(sim_result > superiority_threshold)))

            if (R >= 1e4) { Sys.sleep(time = 60) }

         }

      }

   #-------------------------------------------- Closing up parallel computing --------------------------------------------#

      stopCluster(cl = cl)


#------------------------------------------ Saving image of all the simulations -------------------------------------------#

   save.image('G:/03_Projects/001_Website/jaeger-consulting.github.io/documents/sample_size_allocation_ratio_prior_information.RData')
