path <- "[Directory where files are being saved]"
scriptpath <-  paste0(path, "Scripts/")

temppath <- paste0(path, "Temp/")

source(paste0(scriptpath,"01_selection_simulation_function.R"))


#### Set up parallel structure ####
# Resource:
# https://www.blasbenito.com/post/02_parallelizing_loops_with_r/

library(foreach)
library(doParallel)

library(abind)
acomb <- function(...) abind(..., along=3)

##********************************##
#### SELECT SIMULATION SETTINGS ####
##********************************##

#make simulation settings in the form of a list

N <- 10000 #sample size

inf_base_list <- list(0.7)  #baseline probability of infection

#this time, I'll actually do like the design of the figures, to set the marginal
#probabilities of infection and the confounder, and then induce the difference
#using a risk ratio
agecat_base_list <- list(0.55,0.55)
inf_beta_agecat_list <- list(log(1.7),log(0.65))


pasc_base_list <- list(0.2) #baseline probability of PASC among infected

pasc_beta_agecat_list <- list(0) #log-odds of pasc by agecat among infected

symptom_base_list <- list(0.2) #baseline prevalence of symptom among non-pasc

symptom_betas_pasc_temp <- log(c(rep(1.7,4),rep(1.5,4),rep(1.3,4),rep(1,28))) #log-odds ratio for symptom by pasc status
symptom_betas_pasc_list <- list(symptom_betas_pasc_temp,
                                1.2 * symptom_betas_pasc_temp,
                                1.4 * symptom_betas_pasc_temp
                                )

symptom_betas_agecat_list <- list(numeric(40),
                                  1.7*log(c(rep(1,12),
                                        1.5,1.35,1.20,1.2,
                                        1.5,1.35,1.35,1.2,
                                        1.5,1.50,1.35,1.2,
                                        rep(1,16))), #only on symptoms not associated with PASC
                                  1.7*log(c(1.5,1.35,1.20,1.2,
                                        1.5,1.35,1.35,1.2,
                                        1.5,1.50,1.35,1.2,
                                        rep(1,28)))) #only on symptoms associated with PASC

#I did a whole permutation thing to make sure that the groups were sort of spread out
#with setting 2 being group-sparse, and setting 3 giving each group a mix of zeros and non-zeros
symptom_groups_list <- list(1:40,
                            c(1,3,4,4,2,3,4,4,2,3,4,4,5,5,5,6,6,6,7,7,7,7,7,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9),
                            c(1,2,3,4,5,6,7,8,5,6,7,4,8,1,2,2,3,3,3,3,3,4,5,5,5,6,6,6,6,6,6,6,6,7,7,7,7,7,8,9))

# exp(symptom_betas_agecat_list[[3]] + symptom_betas_pasc_list[[1]]) * symptom_base_list[[1]]

# #make sure that the group sizes are the same across both specifications
# sort(table(symptom_groups_list[[2]]));sort(table(symptom_groups_list[[3]]))
# #check out the associations by symptom by group
# cbind(symptom_groups_list[[2]],symptom_betas_pasc_temp,symptom_betas_agecat_list[[3]])[order(symptom_groups_list[[2]]),]
# cbind(symptom_groups_list[[3]],symptom_betas_pasc_temp,symptom_betas_agecat_list[[3]])[order(symptom_groups_list[[3]]),]
 
symptom_theta_list <- list(5) #amount of correlation of symptoms within a group

#direct relationship between age and the symptoms

# Set up parallel computation ----

#see how many cores we have working!
parallel::detectCores()
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
  # type = "PSOCK"
)
#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
#check if it is registered (optional)
foreach::getDoParRegistered()
#how many workers are available? (optional)
foreach::getDoParWorkers()


# Perform simulation ----

#create a looping structure to work through a bunch of settings...

#might batch out multiple sets of iterations
iter_vec <- 1:500
n_iter <- length(iter_vec)

count_start <- 1
counter <- 1
for(inf_base_ind in seq_along(inf_base_list)){
  inf_base <- inf_base_list[[inf_base_ind]]
  
  for(inf_beta_agecat_ind in seq_along(inf_beta_agecat_list)){
    agecat_base <- agecat_base_list[[inf_beta_agecat_ind]]
    inf_beta_agecat <- inf_beta_agecat_list[[inf_beta_agecat_ind]]
    
    pr_A_Z0_temp <- inf_base / (agecat_base * (exp(inf_beta_agecat)-1) + 1)
    pr_A_Z1_temp <- pr_A_Z0_temp * exp(inf_beta_agecat)
    inf0_age1prop <- (1-pr_A_Z1_temp) * agecat_base / (1-inf_base)
    inf1_age1prop <- pr_A_Z1_temp * agecat_base / inf_base

      for(pasc_base_ind in seq_along(pasc_base_list)){
        pasc_base <- pasc_base_list[[pasc_base_ind]]
        
        for(pasc_beta_agecat_ind in seq_along(pasc_beta_agecat_list)){
          pasc_beta_agecat <- pasc_beta_agecat_list[[pasc_beta_agecat_ind]]
        
          for(symptom_groups_ind in seq_along(symptom_groups_list)){
            symptom_groups <- symptom_groups_list[[symptom_groups_ind]] 
            symptom_theta <- symptom_theta_list[[1]] #hardcode just one level of correlation
            
            for(symptom_base_ind in seq_along(symptom_base_list)){
              symptom_base <- symptom_base_list[[symptom_base_ind]]
              
              for(symptom_betas_agecat_ind in seq_along(symptom_betas_agecat_list)){
                symptom_betas_agecat <- symptom_betas_agecat_list[[symptom_betas_agecat_ind]]
                
                for(symptom_betas_pasc_ind in seq_along(symptom_betas_pasc_list)){
                  symptom_betas_pasc <- symptom_betas_pasc_list[[symptom_betas_pasc_ind]]
                  
                  stopifnot(length(symptom_betas_pasc) == length(symptom_betas_agecat) &
                            length(symptom_betas_pasc) == length(symptom_groups))
                  
                  print(counter); print(Sys.time())
                  # stop("halt!")
                  
                  # if(counter>=count_start){
                  #   # foreach(i = iter_vec, .final = function(x) NULL) %do% {
                  #   foreach(i = iter_vec, .final = function(x) NULL) %dopar% {
                  #     print(paste0("iteration: ", i, " at ", Sys.time()))
                  #     list_out <- selection_simulation_single_function(N=N, #agecat_base=agecat_base,
                  #                                            inf_base=inf_base,
                  #                                            inf1_age1prop=inf1_age1prop, inf0_age1prop=inf0_age1prop,
                  #                                            pasc_base=pasc_base, pasc_beta_agecat=pasc_beta_agecat,
                  #                                            symptom_base=symptom_base,
                  #                                            symptom_groups=symptom_groups, symptom_theta=symptom_theta,
                  #                                            symptom_betas_agecat=symptom_betas_agecat,
                  #                                            symptom_betas_pasc=symptom_betas_pasc,
                  #                                            symptom_pasc_agecat_interact = 0,
                  #                                            seed = 1233 + i, symptom_gen_type="copula",
                  #                                            verbose = FALSE)
                  #     list_out$truth$symptom_betas_agecat_ind <- symptom_betas_agecat_ind
                  #     list_out$truth$symptom_betas_pasc_ind <- symptom_betas_pasc_ind
                  #     list_out$truth$symptom_groups_ind <- symptom_groups_ind
                  # 
                  #     saveRDS(object = list_out, file = paste0(temppath,
                  #                                              "sim_N_",sprintf("%05d", N),
                  #                                              "_setting_",sprintf("%03d", counter),
                  #                                              "_iter_",sprintf("%04d", i),
                  #                                              "_",
                  #                                              ".RDS"))
                  #     NULL
                  #   }
                  # }
                  
                  
                  counter <- counter + 1
                  
                  
                }
              }
            }
          }
        }
      }
    # }
  }
}

#at the end, close the parallel backend
parallel::stopCluster(cl = my.cluster)
