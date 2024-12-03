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

library(dplyr)




# Set up the few key inputs to generating all of the simulation summaries
N <- 10000 #sample size
iter_vec <- 1:500
n_iter <- length(iter_vec)
n_counters <- 54 #vector of total number of simulation settings


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



# Summarize simulation ----

count_start <- 1

# First, run and save simulation summaries on every dataset ----

#loop through counters here too, and figure out then how to best summarize the results across counters
for(analysis in c(#"unweighted",
                  "classweighted",
                  "balanceweighted",
                  # "unweightedadj",
                  "classweightedadj",
                  "balanceweightedadj",
                  NULL)){
  for(counter in count_start:n_counters){
  print(paste0("analysis ", analysis, " counter iteration: ", counter, " at ", Sys.time()))

  # Generate the simulation summaries ----
  ## https://stackoverflow.com/questions/17570806/parallel-for-loop-with-an-array-as-output
  # print(analysis)
    
  foreach(i = iter_vec, .final = function(x) NULL,
          .packages = c("glmnet","ggplot2","pROC","ROCR")) %dopar% {
            
          print(paste0("iteration: ", i, " at ", Sys.time()))
          list_out <- readRDS(file = paste0(temppath,
                                             "sim_N_",sprintf("%05d", N),
                                             "_setting_",sprintf("%03d", counter),
                                             "_iter_",sprintf("%04d", i),
                                             "_",
                                             ".RDS"))
          temp_out <- suppressMessages(simulation_analysis_single_function(list_out, 
                                        list_out$data, analysis, verbose = FALSE))
          saveRDS(object = temp_out, file = paste0(temppath,
                                                   "results_N_",sprintf("%05d", N),
                                                   "_setting_",sprintf("%03d", counter),
                                                   "_iter_",sprintf("%04d", i),
                                                   "_analysis_",analysis,"_",
                                                   ".RDS"))
          NULL
          }
  }
}






# Second, read them in to report ----
sim_settings_list <- coef_mat_list <- 
  selection_mat_list <-
  vector(length=n_counters, mode="list")

#loop through counters here too, and figure out then how to best summarize the results across counters
for(counter in count_start:n_counters){
  print(paste0("counter iteration: ", counter, " at ", Sys.time()))
  
  selection_mat_list[[counter]] <- coef_mat_list[[counter]] <- 
    list(unweighted=NULL, classweighted=NULL,weighted=NULL)
  
  #read just the first one, to get all the "true" specifications for this simulation setting
  sim_settings_list[[counter]] <- readRDS(file = paste0(temppath,
                                                        "sim_N_",sprintf("%05d", N),
                                                        "_setting_",sprintf("%03d", counter),
                                                        "_iter_",sprintf("%04d", 1),
                                                        "_",
                                                        ".RDS"))$truth
  
  # Generate the simulation summaries ----
  ## https://stackoverflow.com/questions/17570806/parallel-for-loop-with-an-array-as-output
  
  for(analysis in c(#"unweighted", 
                    "classweighted","balanceweighted",
                    "classweightedadj","balanceweightedadj",
                    #"adjusted",
                    NULL)){
    print(analysis)
    #    begin_time <- Sys.time()
    selection_mat_list[[counter]][[analysis]] <- 
      foreach(i = iter_vec, .combine='acomb', .multicombine = TRUE) %dopar% {
           temp_out <- readRDS(file = paste0(temppath,
                                              "results_N_",sprintf("%05d", N),
                                              "_setting_",sprintf("%03d", counter),
                                              "_iter_",sprintf("%04d", i),
                                              "_analysis_",analysis,"_",
                                              ".RDS"))
           temp_out
           }

    coef_mat_list[[counter]][[analysis]] <-
      foreach(i = iter_vec, .combine='acomb', .multicombine = TRUE) %dopar% {
      temp_list <- readRDS(file = paste0(temppath,
                                         "sim_N_",sprintf("%05d", N),
                                         "_setting_",sprintf("%03d", counter),
                                         "_iter_",sprintf("%04d", i),
                                         "_",
                                         ".RDS"))
      temp_list[[paste0("coef_mat_",analysis)]]
      }
  }
}

#at the end, close the parallel backend
parallel::stopCluster(cl = my.cluster)

saveRDS(coef_mat_list,  file = paste0(temppath,"coef_mat_list_N_",sprintf("%05d", N),"_",".RDS"))
saveRDS(sim_settings_list,  file = paste0(temppath,"sim_settings_list_N_",sprintf("%05d", N),"_",".RDS"))
# rm(coef_mat_list)

method_names <- c("full", "full_oracle",
                  "lasso_min","lasso_1se",
                  "elastic_min","elastic_1se",
                  "symp_count",
                  NULL)

#make one massively long dataset
out <- vector(mode = "list", length = n_counters)

for(x in 1:n_counters){ #n_counters
  print(paste0("counter: ",x," ",Sys.time()))
  temp_settings_frame <- data.frame(
    N = sim_settings_list[[x]]$N,
    # agecat_base = sim_settings_list[[x]]$agecat_base,
    inf_base = sim_settings_list[[x]]$inf_base,
    inf1_age1prop = sim_settings_list[[x]]$inf1_age1prop,
    inf0_age1prop = sim_settings_list[[x]]$inf0_age1prop,
    inf_beta_agecat = sim_settings_list[[x]]$inf_beta_agecat,
    pasc_base = sim_settings_list[[x]]$pasc_base,
    pasc_beta_agecat = sim_settings_list[[x]]$pasc_beta_agecat,
    n_symptoms = sim_settings_list[[x]]$n_symptoms,
    symptom_base = sim_settings_list[[x]]$symptom_base,
    symptom_groups_ind = sim_settings_list[[x]]$symptom_groups_ind,
    symptom_theta = sim_settings_list[[x]]$symptom_theta,
    symptom_betas_agecat_ind = sim_settings_list[[x]]$symptom_betas_agecat_ind,
    symptom_betas_pasc_ind = sim_settings_list[[x]]$symptom_betas_pasc_ind)

  temp_list <- list()
  for(y in c(#"unweighted", 
    "classweighted","balanceweighted",
    "classweightedadj","balanceweightedadj",
    #"adjusted",
    NULL)){
    for(z in method_names){
      print(paste0(y," ",z))
      temp_list2 <- vector(mode = "list", length = n_iter)
      for(v in iter_vec){
        temp_list2[[v]] <- data.frame(counter=x,
                                       analysis=y, method=z, iter=v,
                                       metric_name = names(selection_mat_list[[x]][[y]][z,,v]),
                                       metric_value = selection_mat_list[[x]][[y]][z,,v])
      }
      temp_list[[paste0(y,"_",z)]] <- do.call(what=rbind, args=temp_list2)
    }
    selection_mat_list[[x]][[y]] <- NULL #free up some space along the way by deleting behind
  }
  
  out[[x]] <- cbind(temp_settings_frame, do.call(what=rbind, args=temp_list))
}
rm(temp_list,temp_list2)
out_frame <- do.call(what=rbind, args=out)
rm(out)
gc()

saveRDS(out_frame,  file = paste0(temppath,"out_frame_N_",sprintf("%05d", N),"_",".RDS"))

