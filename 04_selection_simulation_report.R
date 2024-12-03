path <- "/Users/htr2/Documents/RECOVER/methods/"
scriptpath <-  paste0(path, "Scripts/")

temppath <- paste0(path, "Temp/")

library(dplyr)
library(rsimsum)
library(ggplot2)
library(xtable)
library(tibble)
library(tidyr)

paired_color_palette <- RColorBrewer::brewer.pal(n = 12, name="Paired")[c(10,9,4,3,6,5)]
color_palette <- c(paired_color_palette[1], # "#1F78B4"
                   paired_color_palette[3], # "#33A02C"
                   paired_color_palette[5], # "#E31A1C"
                   "#FFB000")

N <- 10000 #sample size

#get the outframe from the previous script
out_frame_temp <- readRDS(file = paste0(temppath,"out_frame_N_",sprintf("%05d", N),
                                   "_",".RDS"))
#compute the extra "n_selected" variable I forgot to use, back from the true and false positive rates
out_frame <- out_frame_temp %>% bind_rows(
  left_join(
    out_frame_temp %>% filter(metric_name %in% c("tpr")) %>% mutate(tpn = metric_value * 12) %>%
      select(-metric_name,-metric_value),
    out_frame_temp %>% filter(metric_name %in% c("fpr")) %>% mutate(fpn = metric_value * 28) %>%
      select(-metric_name,-metric_value)
  ) %>% mutate(metric_name="n_selected",metric_value=tpn+fpn,.keep="unused")
)
rm(out_frame_temp)
coef_mat_list <- readRDS(file = paste0(temppath,"coef_mat_list_N_",sprintf("%05d", N),"_",".RDS"))

sim_settings_list <- readRDS(file = paste0(temppath,"sim_settings_list_N_",sprintf("%05d", N),"_",".RDS"))
simsettings <- out_frame %>% #group_by(counter, analysis, method, metric_name) %>%
  select(counter,inf1_age1prop,inf0_age1prop, inf_beta_agecat, pasc_base, pasc_beta_agecat,
         symptom_theta,symptom_groups_ind,
         symptom_betas_agecat_ind,symptom_betas_pasc_ind) %>% distinct
rownames(simsettings) <- NULL

##*********************************##
#### OVERALL SUMMARY TABLE (BIG) ####
##*********************************##

summ_out_frame <- out_frame %>%
  group_by(N, counter, inf0_age1prop,inf_base,inf_beta_agecat,
           pasc_base,pasc_beta_agecat,
           n_symptoms,symptom_base,
           symptom_theta,symptom_groups_ind,
           symptom_betas_agecat_ind,
           symptom_betas_pasc_ind, metric_name, method, analysis) %>%
  summarize(metric_mean=mean(metric_value,na.rm=TRUE),
            metric_sd=sd(metric_value, na.rm=TRUE),
            metric_median=median(metric_value,na.rm=TRUE),
            metric_min = min(metric_value,na.rm=TRUE),
            metric_q25 = quantile(metric_value,probs = 0.25, na.rm=TRUE),
            metric_q75 = quantile(metric_value,probs = 0.75, na.rm=TRUE),
            metric_iqr = metric_q75 - metric_q25,
            metric_max = max(metric_value,na.rm=TRUE),
            metric_minnoout = max(metric_q25 - 1.5*metric_iqr,metric_min),
            metric_maxnoout = min(metric_q75 + 1.5*metric_iqr,metric_max),
            count_na=sum(is.na(metric_value))) %>%
  mutate(metric_ll = metric_mean - 1.96*metric_sd, metric_ul = metric_mean + 1.96*metric_sd) %>%
  #restrict to 'medium' setting
  filter(symptom_betas_pasc_ind==2)

summ_out_frame %>% View

#### SIMPLIFIED SIMULATION PLOTS (CONFOUNDING) ####

comparator_vec <- c("lasso_1se",
                    NULL)
summ_out_frame_temp <- summ_out_frame %>% 
  filter(analysis %in% c("classweighted","balanceweighted"), 
         method %in% comparator_vec) %>%
  mutate(y_min = case_when(metric_name=="auc" ~ 0.5,
                           metric_name=="aucpr" ~ 0,
                           metric_name=="tpr" ~ 0,
                           metric_name=="tnr" ~ 0,
                           TRUE ~ NA),
         y_max = case_when(metric_name=="auc" ~ 1,
                           metric_name=="aucpr" ~ 1,
                           metric_name=="tpr" ~ 1,
                           metric_name=="tnr" ~ 1,
                           TRUE ~ NA),
         inf_beta_fact = factor(as.numeric(inf_beta_agecat>1), 
                            levels=c("1","0"),
                            labels=c("Pos. Z/A Assoc.",
                                     "Neg. Z/A Assoc.")),
         symptom_betas_agecat_fact = factor(symptom_betas_agecat_ind, levels=c("1","2","3"),
                            labels=c("No Z/X Assoc.",
                                     "Z/X Non-Over.",
                                     "Z/X Over.")),
         meth_analy_int = interaction(method,analysis),
         meth_fact = factor(meth_analy_int, levels=c(
                                             "lasso_1se.classweighted",
                                             "lasso_1se.balanceweighted",                                             
                                             NULL),
                            labels=c("Lasso (Unadjusted)",
                                     "Lasso (Balancing Weights)",
                                     NULL)),
         metric_fact = factor(metric_name,
                              levels=c("auc","aucpr",
                                       "wilcox_stat",
                                       "wilcox_diff",
                                       "kendall",
                                       "kendall_nonzero",
                                       "tnr","tpr", "n_selected"),
                              labels=c("Area under ROC Curve","Area under Precision-Recall Curve",
                                       "Wilcoxon Rank Sum Test Statistic",
                                       "Hodges–Lehmann Difference in PASC Score Location",
                                       "Est. vs. True Coef. Rank-Correlation",
                                       "Est. vs. True Coef. Rank-Correlation (True Non-Zero)",
                                       "True Negative Rate","True Positive Rate","Selected")),
         symp_inf_int = interaction(symptom_betas_agecat_ind,as.numeric(inf_beta_agecat>1)),
         symp_inf_int2 = factor(symp_inf_int,
                                 levels=c("1.1","1.0",
                                          "2.1","2.0",
                                          "3.1","3.0"),
                                 labels = c("No Z/X Assoc.\nPos. Z/A Assoc.",
                                            "No Z/X Assoc.\nNeg. Z/A Assoc.",
                                            "Z/X Non-Over.\nPos. Z/A Assoc.",
                                            "Z/X Non-Over.\nNeg. Z/A Assoc.",
                                            "Z/X Overlapping\nPos. Z/A Assoc.",
                                            "Z/X Overlapping\nNeg. Z/A Assoc.")))
color_vec <- color_palette[1:length(comparator_vec)]








#make boxplots
#now, output the 2x2 grids of metric plots to pages of a pdf document
cairo_pdf(filename = paste0("/Users/htr2/Documents/RECOVER/methods/Temp/",
                            "sim_boxplots",
                            "_2024-04-22.pdf"),
          onefile=TRUE, width=8,height=5)

comparator_vec <- c("lasso_1se","symp_count")
ggplot(data=out_frame %>% 
         filter(metric_name %in% c("auc","aucpr","wilcox_stat"),
                analysis == "balanceweighted", 
                pasc_beta_agecat==0, 
                symptom_betas_agecat_ind==1, 
                symptom_betas_pasc_ind==2,
                inf0_age1prop==0.75,
                method %in% comparator_vec) %>%
         mutate(
           corr_fact = factor(symptom_groups_ind, levels=c("1","3","2"),
                              labels=c("No Groups","Non-Group Sparse","Group Sparse")),
           meth_fact = factor(method, levels=c("lasso_1se","symp_count"),
             labels=c("Lasso","Symptom Count")),
           metric_fact = factor(metric_name,
                                levels=c("auc","aucpr",
                                         "wilcox_stat",
                                         "kendall",
                                         "kendall_nonzero",
                                         "tnr","tpr","n_selected"),
                                labels=c("Area under ROC Curve","Area under Precision-Recall Curve",
                                         "Wilcoxon Rank Sum Test Statistic",
                                         "Est. vs. True Coef. Rank-Correlation",
                                         "Est. vs. True Coef. Rank-Correlation (True Non-Zero)",
                                         "True Negative Rate","True Positive Rate","Coefs. Selected")),
           y_min = case_when(metric_name=="auc" ~ 0.5,
                             metric_name=="aucpr" ~ 0,
                             metric_name=="tpr" ~ 0,
                             metric_name=="tnr" ~ 0,
                             TRUE ~ NA),
           y_max = case_when(metric_name=="auc" ~ 1,
                             metric_name=="aucpr" ~ 1,
                             metric_name=="tpr" ~ 1,
                             metric_name=="tnr" ~ 1,
                             TRUE ~ NA)),
       aes(x=corr_fact,fill=meth_fact,
           y=metric_value)) + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(facets = ~metric_fact,scales = "free_y") + theme_classic() +
  geom_blank(aes(y=y_min)) + geom_blank(aes(y=y_max)) + 
  xlab("True PASC Association Strength and Symptom Correlation") + ylab("") +
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.text.x = element_text(size=6))

dev.off()









#make boxplots
#now, output the 2x2 grids of metric plots to pages of a pdf document
cairo_pdf(filename = paste0("/Users/htr2/Documents/RECOVER/methods/Temp/",
                            "sim_boxplots_conf",
                            "_2024-05-30.pdf"),
          onefile=TRUE, width=5,height=7)

ggplot(data=summ_out_frame_temp %>% filter(metric_name %in% c("auc","aucpr","wilcox_stat")),
       aes(x=symp_inf_int2,fill=meth_fact,#col=method,
           group=meth_fact,
           lower  = metric_q25,
           upper  = metric_q75,
           middle = metric_median,
           ymin = metric_minnoout,
           ymax = metric_maxnoout)) + 
  # ymin   = metric_q25 - 1.5 * (metric_q75-metric_q25), # optional
  # ymax   = metric_q75 + 1.5 * (metric_q75-metric_q25))) + 
  geom_boxplot(stat = "identity", position=position_dodge(),
               aes(group=interaction(symp_inf_int2, meth_fact)) # optional
  ) + facet_wrap(facets = ~metric_fact,scales = "free_y",ncol = 1) + theme_classic() +
  geom_blank(aes(y=y_min)) + geom_blank(aes(y=y_max)) + ylab("") + 
  xlab("True PASC Association Strength and Symptom Correlation") + 
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.text.x = element_text(size=6))
# theme(legend.position = "bottom",legend.title = element_blank(),
#       legend.text = element_text(size = 16))

dev.off()



summ_out_frame_temp %>% ungroup %>% select(symp_inf_int2, method, meth_fact, metric_fact, 
                                           metric_name, metric_mean, metric_median, metric_sd,metric_q25,metric_q75) %>% 
  filter(metric_name %in% c("tpr","tnr","n_selected","kendall"), method %in% c("lasso_1se")) %>% 
  mutate(label_out = paste0(metric_fact,", ",meth_fact),
         metric_out = ifelse(metric_name=="n_selected",
                             paste0(round(metric_median,3)," (",round(metric_q25,3),"-",round(metric_q75,3),")"),
                             paste0(round(metric_mean,3)," (",round(metric_sd,3),")")),.keep = "unused") %>% 
  pivot_wider(id_cols = c(symp_inf_int2), names_from = label_out,values_from = metric_out) %>% 
  relocate(symp_inf_int2,contains("Selected"),contains("Positive"),contains("Negative")) %>%
  xtable::xtable(x = .) %>% xtable::print.xtable(x=.,include.rownames = FALSE)






summ_out_frame_temp %>% ungroup %>% select(symp_inf_int2, method, meth_fact, metric_fact, 
                                           metric_name, metric_mean, metric_median, metric_sd,metric_q25,metric_q75) %>% 
  filter(metric_name %in% c("tpr","tnr","n_selected","kendall"), 
         method %in% c("lasso_1se")) %>% 
  mutate(label_out = interaction(symp_inf_int2,meth_fact),
         metric_out = ifelse(metric_name=="n_selected",
                             paste0(round(metric_median,3)," (",round(metric_q25,3),"-",round(metric_q75,3),")"),
                             paste0(round(metric_mean,3)," (",round(metric_sd,3),")")),
         .keep = "unused") %>% 
  # select(label_out,meth_fact,metric_fact,metric_out) %>%
  pivot_wider(id_cols = label_out, names_from = metric_fact,values_from = metric_out) %>% 
  select(label_out,`Selected`,`True Positive Rate`,`True Negative Rate`,everything()) %>%
  xtable::xtable(x = .) %>% xtable::print.xtable(x=.,include.rownames = FALSE)

