
######################################################################

            # Find best gamma for each outcome #


# This file does the following:
# 1. function to read in RDS files from tuning
# 2. reading in all RDS files for different outcomes and create tex tables

# IMPORTANT: The data is not available; we have applied for it from a source (the UKBiobank) and interested parties can also apply to the same source

######################################################################


# set paths
if(dir.exists("")){
  mydir <- ""
  setwd(mydir)
  
  my_out_dir <- "" 
  my_figures_out_dir <- "" 
}



############## FUNCTION TO READ IN TUNING DATA ################

# INPUT: 

  # df_string: table identification string used to read in RDS file 
  # outcome_name: outcome name to be used for row names of table (character)

# OUTPUT: 
  # table_res: .tex table showing the number of rejections for different alpha and gamma


find_best_gamma <- function(df_string, outcome_name) {
  
  
  ########## ... platelet crit ################ 
  
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.1 <- readRDS(paste0(mydir, "RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.1"))
  my_outcome_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.1 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.1"))
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.1 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.1"))
  my_outcome_whitenonbritish_unrelatedgamma_matrix_red_f0.1 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedgamma_matrix_red_f0.1"))
  my_outcome_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.1 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.1"))
  
  
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.2 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.2"))
  my_outcome_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.2 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.2"))
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.2 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.2"))
  my_outcome_whitenonbritish_unrelatedgamma_matrix_red_f0.2 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedgamma_matrix_red_f0.2"))
  my_outcome_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.2 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.2"))
  
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.3 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.3"))
  my_outcome_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.3 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.3"))
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.3 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.3"))
  my_outcome_whitenonbritish_unrelatedgamma_matrix_red_f0.3 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedgamma_matrix_red_f0.3"))
  my_outcome_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.3 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.3"))
  
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.4 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.4"))
  my_outcome_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.4 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.4"))
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.4 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.4"))
  my_outcome_whitenonbritish_unrelatedgamma_matrix_red_f0.4 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedgamma_matrix_red_f0.4"))
  my_outcome_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.4 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.4"))
  
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.5 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.5"))
  my_outcome_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.5 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.5"))
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.5 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.5"))
  my_outcome_whitenonbritish_unrelatedgamma_matrix_red_f0.5 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedgamma_matrix_red_f0.5"))
  my_outcome_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.5 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.5"))
  
  
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.6 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.6"))
  my_outcome_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.6 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.6"))
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.6 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.6"))
  my_outcome_whitenonbritish_unrelatedgamma_matrix_red_f0.6 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedgamma_matrix_red_f0.6"))
  my_outcome_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.6 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.6"))
  
  
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.7 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.7"))
  my_outcome_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.7 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.7"))
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.7 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.7"))
  my_outcome_whitenonbritish_unrelatedgamma_matrix_red_f0.7 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedgamma_matrix_red_f0.7"))
  my_outcome_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.7 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.7"))
  
  
  
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.8 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.8"))
  my_outcome_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.8 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.8"))
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.8 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.8"))
  my_outcome_whitenonbritish_unrelatedgamma_matrix_red_f0.8 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedgamma_matrix_red_f0.8"))
  my_outcome_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.8 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.8"))
  
  
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.9 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.9"))
  my_outcome_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.9 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f0.9"))
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.9 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.9"))
  my_outcome_whitenonbritish_unrelatedgamma_matrix_red_f0.9 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedgamma_matrix_red_f0.9"))
  my_outcome_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.9 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedbest_multiplier_matrix_red_f0.9"))
  
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f1 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f1"))
  my_outcome_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f1 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_evals_total_rejections_filter_gamma_red_f1"))
  my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f1 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f1"))
  my_outcome_whitenonbritish_unrelatedgamma_matrix_red_f1 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedgamma_matrix_red_f1"))
  my_outcome_whitenonbritish_unrelatedbest_multiplier_matrix_red_f1 <- readRDS(paste0(mydir,"RDS_tuning_gamma/", df_string, "_whitenonbritish_unrelatedbest_multiplier_matrix_red_f1"))
  
  
  rej_alpha_best_multiplier <- rbind(my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.1, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.2, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.3, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.4, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.5,
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.6,
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.7, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.8, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f0.9, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_red_f1)
  
  rownames(rej_alpha_best_multiplier) <- paste0("Best Multiplier: ", c(paste0(outcome_name, ": Alpha = 0.1"), 
                                                                       paste0(outcome_name, ": Alpha = 0.2"), 
                                                                       paste0(outcome_name, ": Alpha = 0.3"), 
                                                                       paste0(outcome_name, ": Alpha = 0.4"), 
                                                                       paste0(outcome_name, ": Alpha = 0.5"), 
                                                                       paste0(outcome_name, ": Alpha = 0.6"), 
                                                                       paste0(outcome_name, ": Alpha = 0.7"), 
                                                                       paste0(outcome_name, ": Alpha = 0.8"), 
                                                                       paste0(outcome_name, ": Alpha = 0.9"), 
                                                                       paste0(outcome_name, ": Alpha = 1.0")))
  
  
  rej_alpha_unif_multiplier <- rbind(my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.1,
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.2,
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.3, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.4, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.5, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.6, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.7, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.8, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f0.9, 
                                     my_outcome_whitenonbritish_unrelated_fast_total_rejections_filter_gamma_unif_red_f1)

  rownames(rej_alpha_unif_multiplier) <- c(paste0("$\\alpha = 0.", seq(1, 9), "$"), "$\\alpha = 1.0$")

  alpha_rej <- rbind(rej_alpha_unif_multiplier)
  
  
  my_colnames <- c("$\\alpha / 40$", "$\\alpha / 30$", "$\\alpha / 25$", "$\\alpha / 20$","$\\alpha / 15$", 
                   "$\\alpha / 10$", "$\\alpha / 9$" , "$\\alpha / 8$", "$\\alpha / 7$", "$\\alpha / 6$", "$\\alpha / 5$", "$\\alpha / 4$","$\\alpha / 3$", 
                   "$\\alpha / 2$", "$\\alpha$")
  
  colnames(alpha_rej) <- my_colnames
  
  
  stargazer(alpha_rej,
            type = "latex",
            summary=FALSE, rownames = TRUE,
            digits = 2, covariate.labels = c("FDR level $\\alpha$", my_colnames)) -> table_res
  
  tabular_positions <- grep("tabular", table_res)
  table_res <- table_res[tabular_positions[1]:tabular_positions[2]]
  tabular_positions <- grep("tabular", table_res)
  table_res <- table_res[(tabular_positions[1]+1):(tabular_positions[2]-1)]
  table_res <- c(table_res, "\\insertTableNotes")
  table_res <- gsub("\\textbackslash alpha", "\\alpha", table_res, fixed=TRUE)
  table_res <- gsub("\textbackslash alpha", "\\alpha", table_res, fixed=TRUE)
  table_res <- gsub("\\$", "$", table_res, fixed=TRUE)
  
  return(table_res)
  
  
  
}



############## READ TUNING DATA ################


table_res <- find_best_gamma("height", "Height") 
write(table_res, file=paste0(my_out_dir, "whitenonbritish_unrelated_gamma_tune_rejections_height", ".tex"))


table_res <- find_best_gamma("platelet", "Platelet Count") 
write(table_res, file=paste0(my_out_dir, "whitenonbritish_unrelated_gamma_tune_rejections_platelet_count", ".tex"))


table_res <- find_best_gamma("platelet_crit", "Platelet Crit") 
write(table_res, file=paste0(my_out_dir, "whitenonbritish_unrelated_gamma_tune_rejections_platelet_crit", ".tex"))


table_res <- find_best_gamma("platelet_volume", "Platelet Volume") 
write(table_res, file=paste0(my_out_dir, "whitenonbritish_unrelated_gamma_tune_rejections_platelet_volume", ".tex"))


table_res <- find_best_gamma("platelet_width", "Platelet Width") 
write(table_res, file=paste0(my_out_dir, "whitenonbritish_unrelated_gamma_tune_rejections_platelet_width", ".tex"))


table_res <- find_best_gamma("platelet_platelet_volume_platelet_width_platelet_crit", "Platelet Global") 
write(table_res, file=paste0(my_out_dir, "whitenonbritish_unrelated_gamma_tune_rejections_platelet_global", ".tex"))


################ MIN GAMMA PLOTS WNB ###############


height_whitenonbritish_unrelated_min_gamma <- data.frame(readRDS(paste0(mydir, "RDS/", "height_whitenonbritish_unrelated_min_non_inf_gamma")))  %>% mutate(outcome = "height")


global_whitenonbritish_unrelated_min_gamma <- data.frame(readRDS(paste0(mydir, "RDS/", "platelet_platelet_volume_platelet_width_platelet_crit_whitenonbritish_unrelated_min_non_inf_gamma"))) %>% mutate(outcome = "global")

platelet_volume_whitenonbritish_unrelated_min_gamma <- data.frame(readRDS(paste0(mydir, "RDS/", "platelet_volume_whitenonbritish_unrelated_min_non_inf_gamma"))) %>% mutate(outcome = "platelet_volume")

platelet_whitenonbritish_unrelated_min_gamma <- data.frame(readRDS(paste0(mydir, "RDS/", "platelet_whitenonbritish_unrelated_min_non_inf_gamma"))) %>% mutate(outcome = "platelet")

platelet_width_whitenonbritish_unrelated_min_gamma <- data.frame(readRDS(paste0(mydir, "RDS/", "platelet_width_whitenonbritish_unrelated_min_non_inf_gamma")))  %>% mutate(outcome = "platelet width")
platelet_crit_whitenonbritish_unrelated_min_gamma <- data.frame(readRDS(paste0(mydir, "RDS/", "platelet_crit_whitenonbritish_unrelated_min_non_inf_gamma")))  %>% mutate(outcome = "platelet crit")

all_min_gamma <- rbind(height_whitenonbritish_unrelated_min_gamma,  
                       platelet_volume_whitenonbritish_unrelated_min_gamma, 
                       platelet_whitenonbritish_unrelated_min_gamma, 
                       platelet_width_whitenonbritish_unrelated_min_gamma, 
                       platelet_crit_whitenonbritish_unrelated_min_gamma)


all_min_gamma <- all_min_gamma %>% 
  mutate(outcome = ifelse(outcome == "platelet_volume", "platelet volumne", outcome)) %>% 
  mutate(outcome = str_to_sentence(outcome))

resolution_label_strings <- c("single-SNP", "3 kb", "20 kb", "41 kb", "81 kb", "208 kb", "425 kb")


all_min_gamma %>% ggplot(aes(x = Resolution, y = min_gamma, color = outcome)) + geom_line() + 
  labs(x = "Resolution", y = "Smallest finite gamma", color = "Outcome" ) + 
  scale_x_continuous(breaks=seq(0,6,1), labels = resolution_label_strings) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_colour_viridis_d() -> min_gamma


ggsave(paste0(my_figures_out_dir,"wnb_min_gamma_by_outcome.png"), width = 7, height = 6)



