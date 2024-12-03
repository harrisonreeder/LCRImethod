#RECOVER Methods Paper Test Script
library(glmnet)
# library(sparsegl)
# library(gglasso)
library(copula)

library(ggplot2)
library(pROC)
library(ROCR)

#Function that performs simulation study! based on a bunch of inputs----

selection_simulation_single_function <- function(N, #number of individuals to simulate
                                       #agecat_base, #prevalence of demographic covariate overall
                                       inf_base, #prevalence of infection overall
                                       inf1_age1prop, #proportion of infected with age=1
                                       inf0_age1prop, #proportion of uninfected with age=1
                                       pasc_base, #prevalence of pasc among infected
                                       pasc_beta_agecat, #risk ratio between pasc and covariate
                                       n_symptoms, #number of symptoms
                                       symptom_base, #base prevalence of symptom
                                       symptom_groups, #vector "grouping" symptoms
                                       symptom_theta, #copula dependence parameter
                                       symptom_betas_agecat, #risk ratio between symptom and covariate
                                       symptom_betas_pasc, #risk ratio between symptom and pasc
                                       symptom_pasc_agecat_interact, #indicator for whether symptom_betas_agecat are "interaction" coefficients that only affect agecat==1 group
                                       seed = NULL, 
                                       symptom_gen_type="copula",verbose = FALSE){
  # browser()
  
  n_symptoms <- length(symptom_betas_pasc)
  stopifnot(length(symptom_groups) == n_symptoms)
  stopifnot(length(symptom_betas_agecat) == n_symptoms)
  stopifnot(length(symptom_base) %in% c(1, n_symptoms))
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  #vector of names for each method
  method_names <- c(#"full_pasc", "full_oracle_pasc",
    #"lasso_min_pasc","lasso_1se_pasc",
    # "elastic_min_pasc","elastic_1se_pasc",
    # "univariate_pasc",
    
    "full", "full_oracle",  #"forward_aic",
    "lasso_min","lasso_1se",
    "elastic_min","elastic_1se",
    #"univariate", #"symp_count", #this is "done" in the analysis function
    NULL)  
  
  #Generate data----
  ##****************
  
  for(data_iter in c(#"external", 
    "internal")){
    
    #1. generate infection status
    inf_prob <- inf_base
    
    # inf_intercept <- log(inf_base)
    # inf_prob <- exp(inf_intercept + inf_beta_agecat * agecat) #log-linear (risk ratio) scale
    # inf_intercept <- qlogis(inf_base)
    # inf_prob <- plogis(inf_intercept + inf_beta_agecat * agecat)
    
    infection <- rbinom(n=N, size = 1, prob = inf_prob)

    #2. simulate demographic covariate, like age, such that it's more common among infected
    agecat <- rep(NA,N)
    agecat[infection==1] <- rbinom(n=sum(infection==1), size=1, prob=inf1_age1prop)
    agecat[infection==0] <- rbinom(n=sum(infection==0), size=1, prob=inf0_age1prop)

    #3. simulate PASC statuses from a bernoulli (for now, binary PASC vs non-PASC)
    
    pasc_intercept <- log(pasc_base)
    pasc_prob <- exp(pasc_intercept + pasc_beta_agecat * agecat) #log-linear (risk ratio) scale
    # pasc_intercept <- qlogis(pasc_base)
    # pasc_prob <- plogis(pasc_intercept + pasc_beta_agecat * agecat)
    
    pasc_prob[infection==0] <- 0 #uninfected cannot get pasc
    pasc <- rbinom(n=N, size = 1, prob = pasc_prob)
    
    table(infection,pasc)
    
    #4. simulate symptoms based on infection, pasc, and agecat
    
    #vector with number of symptoms in each symptom group
    symptom_group_counts <- as.vector(table(symptom_groups))
    n_groups <- max(symptom_groups)
    symptom_intercepts <- log(symptom_base)
    
    if(symptom_pasc_agecat_interact==1){
      agecat_reg <- agecat * pasc
    } else{
      agecat_reg <- agecat
    }
    
    
    if(symptom_gen_type == "logit"){
      #exclude groups with only one symptom from getting a random effect
      symptom_betas_re <- rep(1,n_symptoms)
      symptom_betas_re[symptom_group_counts[symptom_groups]==1] <- 0
      symptom_individual_re <- MASS::mvrnorm(n = N, mu = numeric(n_groups), 
                                             Sigma = symptom_theta * diag(n_groups))
      symptom_lp_mat <- 
        matrix(data=symptom_intercepts, nrow=N,ncol=n_symptoms, byrow=TRUE) + 
        agecat_reg %*% t(symptom_betas_agecat) + 
        pasc %*% t(symptom_betas_pasc) + 
        t( t(symptom_individual_re[,symptom_groups]) * symptom_betas_re)
      symptom_prob_mat <- plogis(symptom_lp_mat)
      symptom_mat <- sapply(1:n_symptoms, function(j) rbinom(n=N, size=1, prob=symptom_prob_mat[,j]))
    } else{
      #I want this to be an n by p (number of symptoms) matrix of marginal probabilities
      symptom_prob_mat <- exp(
        matrix(data=symptom_intercepts, nrow=N,ncol=n_symptoms, byrow=TRUE) + 
          agecat_reg %*% t(symptom_betas_agecat) + 
          pasc %*% t(symptom_betas_pasc))
      symptom_mat <- do.call(what = cbind, 
         args = lapply(1:n_groups, 
           function(x){
             #step one: generate vector of random variables each on [0,1] under clayton copula
             if(symptom_group_counts[x]>1){
               U <- rCopula(n=N, copula = claytonCopula(param = symptom_theta, 
                                                        dim = symptom_group_counts[x]))
             } else {
               U <- as.matrix(runif(n=N))
             }
             #step two: generate final binary variables using inverse probability transforms
             1 * (U > 1 - symptom_prob_mat[,symptom_groups==x,drop=FALSE])
           }))
    }
    
    if(verbose){ 
      # colMeans(symptom_mat[pasc==0,])
      # colMeans(symptom_mat[pasc==1,])
      
      cbind(colMeans(symptom_mat[pasc==1,])/colMeans(symptom_mat[pasc==0,]),
            exp(symptom_betas_pasc))
      
      cbind(colMeans(symptom_mat[pasc==1,])/(1-colMeans(symptom_mat[pasc==1,]))/
              (colMeans(symptom_mat[pasc==0,])/(1-colMeans(symptom_mat[pasc==0,]))),
            exp(symptom_betas_pasc))
      
      cor_symptom_mat <- cor(symptom_mat,method = "pearson")
      round(cor_symptom_mat,2)
      1 * (abs(cor_symptom_mat)>0.05)
      corrplot::corrplot(corr=cor_symptom_mat)
    }
    
    if(data_iter == "external"){
      temp_data_external <- data.frame(infection,pasc,#age,
                                       agecat,symptom=symptom_mat)
      colnames(temp_data_external) <- sub("\\.",replacement = "",colnames(temp_data_external))
    } else{
      temp_data <- data.frame(infection,pasc,#age,
                              agecat,symptom=symptom_mat)
      colnames(temp_data) <- sub("\\.",replacement = "",colnames(temp_data))
    }
    #second iteration is the one that "persists" below, so that's the data used for fitting everything
    #first iteration is just an extra "external" dataset to be used for assessment of external predictive performance
  }
  
  #Run analyses ----
  ##****************
  
  temp_list <- list()
  
  #class imbalance weight overall
  non_weights <- rep(1,N)
  overall_weights <- rep(1,N)
  overall_weights[infection==0] <- sum(infection==1) / sum(infection==0)
  #class imbalance weight stratified by demographic covariate
  agecat_weights <- rep(1,N)
  agecat_weights[infection==0 & agecat==1] <- sum(infection[agecat==1]==1) / sum(infection[agecat==1]==0)
  agecat_weights[infection==0 & agecat==0] <- sum(infection[agecat==0]==1) / sum(infection[agecat==0]==0)
  
  table(agecat,agecat_weights,infection)
  #there are four sets of analyses: 
  # "unweighted" just using the method directly (without class imbalance weighting), 
  # "classweighted" using class imbalance weighting,
  # "balanceweighted" using a balancing weight defined on agecat, and
  # "both" using a balancing weight defined on agecat, and adjusting for the covariate directly
  
  # analysis <- "adjusted"
  for(analysis in c(#"unweighted", 
                    "classweighted", 
                    "balanceweighted",
                    #"unweightedadj",
                    "classweightedadj",
                    "balanceweightedadj",
                    NULL)){
    
    null_form <- infection ~ 1
    
    full_form <- as.formula(paste0("infection ~ ",paste(paste0("symptom",1:n_symptoms),collapse = " + ")))
    if(analysis %in% c("unweightedadj","classweightedadj","balanceweightedadj")){ full_form <- update(full_form, . ~ . + agecat)}
    
    oracle_form <- as.formula(paste0("infection ~ ",paste(paste0("symptom",which(symptom_betas_pasc != 0)),collapse = " + ")))
    if(analysis %in% c("unweightedadj","classweightedadj","balanceweightedadj")){ oracle_form <- update(oracle_form, . ~ . + agecat)}
    
    if(analysis %in% c("unweighted","unweightedadj")){
      weights_temp <- non_weights
    } else if(analysis %in% c("classweighted","classweightedadj")){
      weights_temp <- overall_weights
    } else if(analysis %in% c("balanceweighted","balanceweightedadj")){
      weights_temp <- agecat_weights
    } else { 
      stop("something went wrong with the analysis specifications")
    }
    
    null_fit <- glm(formula = null_form, data=temp_data, weights = weights_temp, 
                    family = binomial(link="logit"))
    full_fit <- glm(formula = full_form, data=temp_data, weights = weights_temp, 
                    family = binomial(link="logit"), model=FALSE,x=FALSE,y=FALSE)
    if(verbose) print(summary(full_fit))
    
    full_oracle_fit <- glm(formula = oracle_form, data=temp_data, weights = weights_temp, 
                           family = binomial(link="logit"), model=FALSE,x=FALSE,y=FALSE)
    if(verbose) print(summary(oracle_fit))
    # round(vcov(full_fit),5)

    ##now, to do all of the penalized models##
    ##**************************************##
    
    temp_symptom_mat <- if(analysis %in% c("unweightedadj", "classweightedadj","balanceweightedadj")) cbind(symptom_mat,agecat) else symptom_mat
    colnames(temp_symptom_mat) <- c(paste0("symptom",1:n_symptoms),  
        if(analysis %in% c("unweightedadj", "classweightedadj","balanceweightedadj")) "agecat") 
    
    cv_lasso_fit <- cv.glmnet(x = temp_symptom_mat, y = infection, family="binomial", 
                              weights = weights_temp, 
                              type.measure = "class", nfolds=10, alpha=1)
    if(verbose){
      print(cv_lasso_fit)
      plot(cv_lasso_fit,main = paste0("lasso_",analysis))
      print(coef(cv_lasso_fit, s = "lambda.min"))
      print(coef(cv_lasso_fit, s = "lambda.1se"))
    }
    
    cv_elastic_fit <- cv.glmnet(x = temp_symptom_mat, y = infection, family="binomial", 
                                weights = weights_temp, 
                                type.measure = "class", nfolds=10, alpha=0.5)
    
    if(verbose){
      print(cv_elastic_fit)
      plot(cv_elastic_fit,main = paste0("lasso_",analysis))
      print(coef(cv_elastic_fit, s = "lambda.min"))
      print(coef(cv_elastic_fit, s = "lambda.1se"))
    }
    
    # cv_lasso_pasc_fit <- cv.glmnet(x = temp_symptom_mat, y = pasc, family="binomial", weights = weights_temp_pasc, 
    #                           type.measure = "class", nfolds=10, alpha=1)
    # if(verbose){
    #   print(cv_lasso_pasc_fit)
    #   plot(cv_lasso_pasc_fit,main = paste0("lasso_",analysis))
    #   print(coef(cv_lasso_pasc_fit, s = "lambda.min"))
    #   print(coef(cv_lasso_pasc_fit, s = "lambda.1se"))
    # }
    # 
    # cv_elastic_pasc_fit <- cv.glmnet(x = temp_symptom_mat, y = pasc, family="binomial", weights = weights_temp_pasc, 
    #                                type.measure = "class", nfolds=10, alpha=0.5)
    # if(verbose){
    #   print(cv_elastic_pasc_fit)
    #   plot(cv_elastic_pasc_fit,main = paste0("lasso_",analysis))
    #   print(coef(cv_elastic_pasc_fit, s = "lambda.min"))
    #   print(coef(cv_elastic_pasc_fit, s = "lambda.1se"))
    # }
    
    #Extract coefficients ----
    ##*************************************************
    
    coef_full <- coef(full_fit)
    
    #slot the "oracle" covariates into a full covariate vector, with 0s elsewhere
    coef_oracle_temp <- coef(full_oracle_fit)
    coef_oracle_temp2 <- numeric(length(coef_full))
    names(coef_oracle_temp2) <- names(coef_full)
    coef_oracle_temp2[names(coef_oracle_temp)] <- coef_oracle_temp
    
    # #slot the pasc "oracle" covariates into a full covariate vector, with 0s elsewhere
    # coef_oracle_pasc_temp <- coef(full_oracle_pasc_fit)
    # coef_oracle_pasc_temp2 <- numeric(length(coef_full))
    # names(coef_oracle_pasc_temp2) <- names(coef_full)
    # coef_oracle_pasc_temp2[names(coef_oracle_pasc_temp)] <- coef_oracle_pasc_temp
    
    #matrix of coefficients for each symptom
    coef_mat <- rbind(
      cbind(
        #full_pasc=coef(full_pasc_fit),
        #full_oracle_pasc = coef_oracle_pasc_temp2,
        #lasso_min_pasc=as.numeric(coef(cv_lasso_pasc_fit, s = "lambda.min")),
        #lasso_1se_pasc=as.numeric(coef(cv_lasso_pasc_fit, s = "lambda.1se")),
        #elastic_min_pasc=as.numeric(coef(cv_lasso_pasc_fit, s = "lambda.min")),
        #elastic_1se_pasc=as.numeric(coef(cv_lasso_pasc_fit, s = "lambda.1se")),
        # univariate_pasc = c("(Intercept)"=NA, univariate_pasc_fits[,1] * as.numeric(univariate_pasc_fits[,2] < 0.05),
        #                if(analysis %in% c("adjusted", "both")) NA),
        
        full=coef_full,
        full_oracle=coef_oracle_temp2,
        lasso_min=as.numeric(coef(cv_lasso_fit, s = "lambda.min")),
        lasso_1se=as.numeric(coef(cv_lasso_fit, s = "lambda.1se")),
        elastic_min=as.numeric(coef(cv_elastic_fit, s = "lambda.min")),
        elastic_1se=as.numeric(coef(cv_elastic_fit, s = "lambda.1se")),
        # univariate = c("(Intercept)"=NA, univariate_fits[,1] * as.numeric(univariate_fits[,2] < 0.05),
        #                if(analysis %in% c("adjusted", "both")) NA),
        NULL),
      agecat = if(analysis %in% c("unweighted","classweighted","balanceweighted")) 0 else NULL) #make it so that "agecat" coefficient always exists
    
    stopifnot(colnames(coef_mat) == method_names)
    
    temp_list[["truth"]] <- 
      list(N=N, #agecat_base=agecat_base,
           inf1_age1prop=inf1_age1prop,inf0_age1prop=inf0_age1prop,
           inf_beta_agecat = inf1_age1prop / inf0_age1prop, #risk ratio implied by the raw proportions
           inf_base=inf_base,
           pasc_base=pasc_base, pasc_beta_agecat=pasc_beta_agecat,
           n_symptoms=n_symptoms, symptom_base=symptom_base,
           symptom_groups=symptom_groups, symptom_theta = symptom_theta,
           symptom_betas_agecat=symptom_betas_agecat, 
           symptom_betas_pasc=symptom_betas_pasc,
           symptom_pasc_agecat_interact=symptom_pasc_agecat_interact) #fill this with all of the true values of the various parameters
    temp_list[["data"]] <- temp_data
    temp_list[[paste0("coef_mat_",analysis)]] <- coef_mat
  }
  
  temp_list
  
}



#Function to "analyze" a single simulation ----

simulation_analysis_single_function <- function(temp_list, temp_data, analysis, verbose=TRUE){
  
  # browser()
  
  #vector of names for each method
  method_names <- c(# "full_pasc", "full_oracle_pasc",
    # "lasso_min_pasc","lasso_1se_pasc",
    # "elastic_min_pasc","elastic_1se_pasc",
    # "univariate_pasc",
    
    "full", "full_oracle",
    "lasso_min","lasso_1se",
    "elastic_min","elastic_1se",
    # "univariate",
    "symp_count",
    NULL)
  
  n_symptoms <- temp_list$truth$n_symptoms
  symptom_betas_pasc <- temp_list$truth$symptom_betas_pasc
  
  #make a coef_mat with just the symptom coefficients, no covariates
  coef_mat_symptoms <- 
    cbind(temp_list[[paste0("coef_mat_",analysis)]][paste0("symptom",1:n_symptoms),,drop=FALSE],
          symp_count = rep(1,n_symptoms))
  
  #Generate predictions ----
  ##************************
  
  #Get the final symptom matrix (leaving off agecat)
  symptom_mat_noagecat <- as.matrix(temp_data[,paste0("symptom",1:n_symptoms),drop=FALSE])
  
  #now, create a matrix of "PASC Scores" as defined by JAMA paper
  #adding the "dumb" comparator of just counting number of symptoms
  #following JAMA paper, round to one decimal when deciding selection
  pred_noint_mat <- cbind(symptom_mat_noagecat %*% coef_mat_symptoms)
  
  # #rescale columns to be between 0 and 10, just for simplicity
  # pred_noint_mat_rescaled <- apply(pred_noint_mat, MARGIN = 2, 
  #                                  FUN = scales::rescale, to = c(0,10))
  
  #rescale columns to have mean 0 and sd 1
  # pred_noint_mat_rescaled <- scale(pred_noint_mat)
  
  #if all coefficients are 0, still scale everything else!
  pred_noint_mat_rescaled <- apply(pred_noint_mat, MARGIN=2,
                                   FUN = function(x){
                                     mean_temp <- mean(x)
                                     centered <- x - mean_temp
                                     sd_temp <- sd(centered)
                                     if(sd_temp != 0){
                                       return(centered/sd_temp)
                                     } else{
                                       return(centered)
                                     }
                                   })
  
  # Plot score distributions ----
  ##*****************************
  
  if(verbose){
    #convert predictions into a long-format dataset for plotting
    pred_frame_temp <- data.frame(id=rep(1:N,length(method_names)),
                                  fit_type=factor(rep(method_names,each=N),levels = method_names),
                                  pred_noint=as.vector(pred_noint_mat),
                                  pred_noint_rescaled=as.vector(pred_noint_mat_rescaled))
    pred_frame <- cbind(pred_frame_temp,
                        agecat=temp_data$agecat,
                        infection=temp_data$infection,
                        pasc=temp_data$pasc, analysis=analysis)
    
    #boxplots of score distribution by model, stratified by true pasc
    print(ggplot(data=pred_frame,mapping=aes(x=pred_noint, y=fit_type, fill=as.factor(pasc))) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(height=0.4,width=0,alpha=0.005))
    print(ggplot(data=pred_frame,mapping=aes(x=pred_noint_rescaled, y=fit_type, fill=as.factor(pasc))) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(height=0.4,width=0,alpha=0.005))
    # #actual predicted probability
    # ggplot(data=pred_frame,mapping=aes(x=expit(pred), y=fit_type, fill=as.factor(pasc))) +
    #   geom_boxplot(outlier.shape = NA) +
    #   geom_jitter(height=0.4,width=0,alpha=0.02)
    
    # #they're really just like nearly identical that's the issue...
    #attempt to see if a ridgeplot would be more useful, but the answer is basically no
    # ggplot(data=pred_frame,mapping=aes(x=pred_noint_rescaled, after_stat(density), fill=as.factor(pasc))) +
    #   geom_histogram(position="identity",alpha=0.35) + facet_grid(rows=vars(fit_type)) + theme_classic()
    # print(ggplot(data=pred_frame,mapping=aes(x=pred_noint_rescaled, y=fit_type, fill=as.factor(pasc))) +
    #   geom_density_ridges(stat = "binline", bins=30, position = "identity", alpha=0.75))
    
  }
  
  # PASC Score Performance Metrics ----
  ##***********************************
  
  # estimate OR of true PASC vs. PASC Score (rescaled) ----
  ##*******************************************************
  
  score_assoc_mat <- matrix(nrow=length(method_names), ncol = 10, 
                            dimnames = list(method_names, c("logor","logor_se","logor_p","logor_ll","logor_ul",
                                                            "wilcox_stat","wilcox_diff","wilcox_ll","wilcox_ul","wilcox_p")))
  for(method in method_names){
    #rescale all the scores so that they range from 0 to 10 (though they keep their internal distribution!)
    # temp <- pred_noint_mat[,method]
    score_temp <- pred_noint_mat_rescaled[,method]
    pasc_temp <- temp_data$pasc
    
    temp_fit <- glm(pasc_temp ~ score_temp, family=binomial(link="logit"), singular.ok = TRUE)
    temp_coef <- tryCatch(expr = coef(summary(temp_fit))["score_temp",c("Estimate","Std. Error","Pr(>|z|)")],
                          error = function(e) return(c(NA,NA,NA)))
    temp_ci <- temp_coef[1] + qt(p = c(0.025,0.975),df = N) * temp_coef[2]
    
    suppressWarnings(temp_wilcox <- tryCatch(expr = wilcox.test(score_temp[pasc_temp==1], score_temp[pasc_temp==0], paired = FALSE, conf.int=TRUE),
                            error = function(e) return(NULL)))
    temp_wilcox <- if(!is.null(temp_wilcox)) c(temp_wilcox$statistic,temp_wilcox$estimate,temp_wilcox$conf.int,temp_wilcox$p.value) else rep(NA,5)
    score_assoc_mat[method,] <- c(temp_coef,temp_ci,temp_wilcox)
    
  }
  
  #flag which did and didn't fit successfully
  method_names_temp_lgl <- !is.na(score_assoc_mat[,"logor"])
  
  if(verbose){
    #finally, a few plots showing the relationship between (rescaled) PASC score and true pasc probability
    #these are logistic regression and logistic GAM models, plotted with logit y to show
    #the differences between the models more clearly
    print(ggplot(data = pred_frame, mapping=aes(x=pred_noint_rescaled, y = pasc, color=fit_type)) +
            geom_smooth(se = FALSE, method = "glm",  method.args = list(family = "binomial")) +
            # geom_point(data=pred_frame[pred_frame$fit_type=="full",],
            #            mapping=aes(x=pred_noint_rescaled, y = pasc), alpha=0.02) +
            theme_classic() + coord_trans(y="logit") +
            scale_color_manual(values = color_palette[which(method_names_temp_lgl)]))
    print(ggplot(data = pred_frame, mapping=aes(x=pred_noint_rescaled, y = pasc, color=fit_type)) +
            geom_smooth(se = FALSE, method = "gam",  method.args = list(family = "binomial")) +
            theme_classic() + coord_trans(y="logit") +
            scale_color_manual(values = color_palette[which(method_names_temp_lgl)]))
    #if I wanted to save the underlying data from these plots I could...
    # #https://stackoverflow.com/questions/9789871/method-to-extract-stat-smooth-line-fit
  }
  
  # estimated ROC curve and AUC-ROC values, with 95% CI estimate ----
  ##*****************************************************************
  
  #lots of options given here: https://rviews.rstudio.com/2019/03/01/some-r-packages-for-roc-curves/
  
  roc_out <- lapply(method_names, function(i) pROC::roc(response=temp_data$pasc, predictor=pred_noint_mat[,i], ci=TRUE))
  auc_out <- t(sapply(roc_out,function(x) x$ci[c(2,1,3)]))
  dimnames(auc_out) <- list(method_names,c("auc","auc_ll","auc_ul"))
  auc_out
  
  #estimate AUC-PR values ----
  ##**************************
  
  roc_out2 <- lapply(method_names, function(i) prediction(pred_noint_mat[,i], temp_data$pasc))
  aucpr_out <- sapply(roc_out2, function(x) performance(x,measure="aucpr")@y.values[[1]])
  names(aucpr_out) <- method_names
  aucpr_out
  
  if(verbose){
    plot(performance(roc_out2[[1]], measure = "prec", x.measure = "rec"), col=color_palette[1], lwd=2)
    for(method_ind in 2:length(method_names)){
      plot(performance(roc_out2[[method_ind]], measure = "prec", x.measure = "rec"),
           add=TRUE,col=color_palette[method_ind], lwd=2)
    }
    legend(x="topright",legend=method_names,fill=color_palette[1:length(method_names)],cex=0.8)
    
    #plot PR curve
    plot(performance(roc_out2[[1]], measure = "tpr", x.measure = "fpr"), col=color_palette[1], lwd=2)
    for(method_ind in 2:length(method_names)){
      plot(performance(roc_out2[[method_ind]], measure = "tpr", x.measure = "fpr"),
           add=TRUE,col=color_palette[method_ind], lwd=2)
    }
    legend(x="bottomright",legend=method_names,fill=color_palette[1:length(method_names)],cex=0.8)
    
    # plot(performance(roc_out2[[1]], measure = "lift", x.measure = "rpp"), col=color_palette[1], lwd=2)
    # for(method_ind in 2:length(method_names)){
    #   plot(performance(roc_out2[[method_ind]], measure = "lift", x.measure = "rpp"),
    #        add=TRUE,col=color_palette[method_ind], lwd=2)
    # }
    # legend(x="topright",legend=method_names,fill=color_palette[1:length(method_names)],cex=0.8)
  }
  
  pred_performance_mat <- cbind(auc_out, #column names already included
                                aucpr = aucpr_out,
                                score_assoc_mat)
  
  #Compute and save selection metrics ----
  ##*****************************
  
  # #precompute rank order metrics
  # #matrix giving the rank order from largest to smallest estimates, by column
  # rank_mat <- apply(X = -coef_mat_symptoms, MARGIN = 2, rank, ties.method = "max")
  # #comparing the rank order to the "true" rank order, where if there is a true tie, then as long as 
  # #the estimates are correctly situated relative to the other estimates, then within the tied group order doesn't matter
  # #we also reshuffle rows so the truly largest symptoms are at the top
  # rank_mat_compared <- (rank_mat <= rank(-symptom_betas_pasc, ties.method = "max"))[order(symptom_betas_pasc, decreasing = TRUE),]
  # #next, we sum down each column, to see how many symptoms are in the correct order
  # rank_mat_compared2 <- apply(rank_mat_compared, 2, cumsum)
  # #finally, look from the top to count how many symptoms in a row are correctly ordered and give that number
  # rank_mat_compared3 <- apply(rank_mat_compared2, 2, function(x) x == seq_along(x))
  # #as a buffer, how about saying the largest x such that your top ranked x contain at least x-1 of the top symptoms
  # rank_mat_compared4 <- apply(rank_mat_compared2, 2, function(x) x >= seq_along(x)-1)
  
  #visually assess the sorted matrix that the rank thing is actually representing
  if(verbose) print(cbind(true=symptom_betas_pasc,coef_mat_symptoms)[order(symptom_betas_pasc,decreasing=TRUE),])
  
  nonzero_ind <- which(symptom_betas_pasc !=0)
  
  #finally, put all of the metrics together!
  suppressWarnings(
    selection_mat <- rbind(
      cbind(
        tpr = colSums((coef_mat_symptoms != 0) & (symptom_betas_pasc != 0)) / sum(symptom_betas_pasc != 0),
        tnr = colSums((coef_mat_symptoms == 0) & (symptom_betas_pasc == 0)) / sum(symptom_betas_pasc == 0),
        fpr = colSums((coef_mat_symptoms != 0) & (symptom_betas_pasc == 0)) / sum(symptom_betas_pasc == 0),
        fnr = colSums((coef_mat_symptoms == 0) & (symptom_betas_pasc != 0)) / sum(symptom_betas_pasc != 0),
        sir = colMeans(sign(coef_mat_symptoms) != sign(symptom_betas_pasc)),
        spearman = apply(X = coef_mat_symptoms, MARGIN = 2, 
                         FUN = function(x) {
                           tryCatch(expr = cor(x = x, y = symptom_betas_pasc, method = "spearman"),
                                    error = function(e) return(NA))}
                         ),
        kendall = apply(X = coef_mat_symptoms, MARGIN = 2, 
                        FUN = function(x) {
                          tryCatch(expr = cor(x = x, y = symptom_betas_pasc, method = "kendall"),
                                   error = function(e) return(NA))}),
        spearman_nonzero = apply(X = coef_mat_symptoms[nonzero_ind,,drop=FALSE], MARGIN = 2, 
                         FUN = function(x) {
                           tryCatch(expr = cor(x = x, y = symptom_betas_pasc[nonzero_ind], method = "spearman"),
                                    error = function(e) return(NA))}
        ),
        kendall_nonzero = apply(X = coef_mat_symptoms[nonzero_ind,,drop=FALSE], MARGIN = 2, 
                        FUN = function(x) {
                          tryCatch(expr = cor(x = x, y = symptom_betas_pasc[nonzero_ind], method = "kendall"),
                                   error = function(e) return(NA))})#,
        # top = colSums(rank_mat_compared3),
        # topm1 = colSums(rank_mat_compared4)
      )
    )
  )
  
  cbind(selection_mat,pred_performance_mat)
  
}
