#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)

######## Read in Arguments ###########

# IMPORTANT: The data is not available; we have applied for it from a source (the UKBiobank) and interested parties can also apply to the same source


library(Surrogate)
# Population here should be the population used for training

population <- as.character(args[1])
if(is.na(population)) population <- "whitenonbritish_unrelated"

outcomes <- as.character(args[2])
outcomes <- eval(parse(text = outcomes)) 
if(is.na(outcomes[1])) outcomes <- "platelet"

alpha <- as.numeric(args[3])
if(is.na(alpha)) alpha <- 0.1

run_from_scratch <- as.character(args[4])
run_from_scratch <- eval(parse(text = run_from_scratch)) 
if(is.na(run_from_scratch)) run_from_scratch <- TRUE

match_matteo_group <- as.character(args[5])
match_matteo_group <- eval(parse(text = match_matteo_group)) 
if(is.na(match_matteo_group)) match_matteo_group <- FALSE

print("This is the population used in this script:")
print(population)

print("This is the outcome used in this script:")
print(outcomes)

print("This is alpha used in this script:")
print(alpha)

print("Is this script running from scratch?")
print(run_from_scratch)

print("Are we trying to match Matteo's lasso groups?")
print(match_matteo_group)

######## Define Defaults, Set Working Directory and Read in Functions ###########
set.seed(2022)


# set paths
if(dir.exists("")){
  mydir <- ""
  setwd(mydir)
  
  my_out_dir <- "" 
}

source(paste0(mydir, "utils/load_aux.R"))


# define a default sequence of gammas to evaluate
gamma_matrix <- matrix(NA, nrow = 7, ncol = 15)

gamma_matrix[, 1] <- c(rep(alpha / 40, 7))
gamma_matrix[, 2] <- c(rep(alpha / 30, 7))
gamma_matrix[, 3] <- c(rep(alpha / 25 , 7))
gamma_matrix[, 4] <- c(rep(alpha / 20, 7))
gamma_matrix[, 5] <- c(rep(alpha / 15, 7))
gamma_matrix[, 6] <- c(rep(alpha / 10, 7))
gamma_matrix[, 7] <- c(rep(alpha / 9, 7))
gamma_matrix[, 8] <- c(rep(alpha / 8, 7))
gamma_matrix[, 9] <- c(rep(alpha / 7, 7))
gamma_matrix[, 10] <- c(rep(alpha / 6, 7))
gamma_matrix[, 11] <- c(rep(alpha / 5, 7))
gamma_matrix[, 12] <- c(rep(alpha / 4, 7))
gamma_matrix[, 13] <- c(rep(alpha / 3, 7))
gamma_matrix[, 14] <- c(rep(alpha / 2, 7))
gamma_matrix[, 15] <- c(rep(alpha , 7))


outcome_string <- paste0(outcomes, collapse = "_")

########## Set up matrices to hold results ######## 
all_resolutions <- seq(0, 6)
chr <- seq(1, 22)

# define matrices to hold the results: 
total_rejections_filter_gamma_frac <- matrix(NA, nrow = 1, ncol = ncol(gamma_matrix)) # hold results for max rejection with multiplier
total_rejections_filter_gamma_eval <- matrix(NA, nrow = 1, ncol = ncol(gamma_matrix)) # hold results for uniform multiplier only 
total_rejections_filter_gamma_unif <- matrix(NA, nrow = 1, ncol = ncol(gamma_matrix)) # hold results for original e-value definition only

best_multiplier_matrix <- matrix(NA, nrow = 7, ncol = ncol(gamma_matrix))
counter <- 0

########### Start loop over all gamma #########

for(gm in seq(1, ncol(gamma_matrix))) {
  
  print("Current gamma column is:")
  print(gamma_matrix[1, gm])
  counter <- counter + 1
  
  # list containing all possible groups, with their fractions, evalues etc
  cand_groups <- list()
  
  # containing snps
  snps <- c()
  
  # containing rejections by global null
  rej_global_res <- c()
  
  # containing rejections by the knockoff
  rej_knockoff <- c()
  
  # information for every snp-group-res combination
  snp_group_chr_matrix <- c()
  
  # smallest gamma that is not infinity
  min_non_inf_gamma <- c()
  
  ################## ...Read in Lasso Results #############
  
  # see if previous files exist
  if(!(run_from_scratch) & file.exists(paste0(mydir, outcome_string, "_", population, "_cand_groups_all_res_gamma_matrix_f", alpha, "g", gamma_matrix[1, gm]))) {
    
    cand_groups <-  readRDS(paste0(mydir, outcome_string, "_", population, "_cand_groups_all_res_gamma_matrix_f", alpha, "g", gamma_matrix[1, gm]))
    snp_group_chr_matrix <- readRDS(paste0(mydir, outcome_string, "_", population, "_snp_group_chr_matrix_gamma_matrix_f", alpha, "g", gamma_matrix[1, gm]))
    snps <- readRDS(paste0(mydir, outcome_string, "_", population, "_snps_gammas_matrix_f", alpha, "g", gamma_matrix[1, gm]))
    rej_global_res <- readRDS(paste0(mydir, outcome_string, "_", population, "_rej_global_res_gamma_matrix_f",alpha, "g", gamma_matrix[1, gm]))
    
  } else {
    
    for(r in all_resolutions) {
      
      print(r)
      
      res <- read_process_UKB_lasso(resolution = r, 
                                     population = population,
                                     alpha = alpha, 
                                     gamma = gamma_matrix[, gm],
                                     mydir = mydir, 
                                     dirpath_group = "UKB_data_input/group_information/ukb_gen_chr",
                                     dirpath_lasso = "UKB_data_input/ukb_stats/",
                                     outdir_path = paste0(my_out_dir, "snp_group_w/"), 
                                     outcomes = outcomes, 
                                     chr = chr, 
                                     smallest_possible_gamma = FALSE, 
                                     gamma_resolution_adaptive = TRUE, 
                                     match_matteo_group = match_matteo_group)
      
      cand_groups <- c(cand_groups, res$cand_groups_resolution)
      snps <- c(snps, res$snps)
      rej_global_res <- rbind(rej_global_res, data.frame(res$rejected_matrix_ebh_global_res))
      snp_group_chr_matrix <- rbind(snp_group_chr_matrix, res$snp_group_chr_matrix_evals)
      rej_knockoff <- rbind(rej_knockoff, res$rejected_knockoff)
      min_non_inf_gamma <- rbind(min_non_inf_gamma, data.frame(min_gamma = res$min_non_inf_gamma) %>% mutate(Resolution = r))
      
      
    }
    
    
    saveRDS(min_non_inf_gamma, file=paste0(my_out_dir, "RDS/",outcome_string, "_", population, "_min_non_inf_gamma")) # this should be the same for every gamma and alpha

    
  }  
  

  ################## ... Data Wrangling of Lasso Results #############
  
  # create chr-group-resolution identifier
  snp_group_chr_matrix$Res_CHR_Group <- paste0(snp_group_chr_matrix$Resolution, "_", snp_group_chr_matrix$CHR_Group)

  #  group size for each resolution-CHR-group
  counts_per_res_chr_group <- snp_group_chr_matrix %>% 
    group_by(Res_CHR_Group) %>%
    count()

  # merge group size into dataset
  snp_group_chr_matrix_wcounts <- dplyr::left_join(snp_group_chr_matrix, counts_per_res_chr_group, 
                                                   by = c("Res_CHR_Group"))
  
  # calculate weight-adjusted fraction
  snp_group_chr_matrix_wcounts <- snp_group_chr_matrix_wcounts %>% 
    mutate(frac_adjusted = fracs / n)
  
  # filter fractions > 0 just to do the adjustment
  snp_group_chr_matrix_wcounts_filtered_nonzero <- snp_group_chr_matrix_wcounts %>% filter(fracs > 0)
  
  # those are the ones with zero fraction and zero evalue
  zero_fracs <- snp_group_chr_matrix_wcounts %>% filter(fracs == 0) 
  
  ###### ... Pre-Filter Results to best in block######
  
  if(!(run_from_scratch) & file.exists(paste0(mydir, outcome_string, "_", population, "_prev_gamma_matrix_f", alpha, "g", gamma_matrix[1, gm]))) {
    
    prev <- readRDS(paste0(mydir, outcome_string, "_", population, "_prev_gamma_matrix_f", alpha, "g", gamma_matrix[1, gm]))
    
  } else {
    
  
      # COMPARE RESOLUTION 1 TO RESOLUTION 0 USING FILTERED SET WITH NONZERO FRACTION
      compare_0_vs_1 <- return_selected_for_two_res(snp_group_chr_matrix_wcounts_filtered_nonzero, bigger_resolution = 1, smaller_resolution = 0)
      prev <- compare_0_vs_1
      
      # STEPWISE COMPARISON USING FILTERED SET WITH NONZERO FRACTION
      # Get the "kept" groups from the previous step, add those from the next resolution
      # run the comparison as described above
      for(r in 2:6) {
        
        
        print(r)
        next_resolution <- snp_group_chr_matrix_wcounts_filtered_nonzero %>% filter(Resolution == r)
        combined_prev_with_next <- rbind(prev, next_resolution)
        
        prev <- return_selected_for_two_res(combined_prev_with_next, bigger_resolution = r, smaller_resolution = r-1)
        
      }

  }
  
  # is every snp included only once? (Check that location worked)
  length(prev$SNP) == length(unique(prev$SNP))
  
  # rename to avoid confusion
  filtered_fracs <- prev 
  
  # total number groups: based on entire set 
  total_number_groups <- length(unique(snp_group_chr_matrix  %>% pull(Res_CHR_Group)))
  # number filtered groups: based on filtered set
  number_filtered_groups <- length(unique(filtered_fracs  %>% pull(Res_CHR_Group)))
  # number groups with zero fraction
  total_number_groups_zero_frac <- length(unique(zero_fracs$Res_CHR_Group))
  
  # the number of groups with nonzero fraction who were filtered out + the number of selected groups + zero groups = total number groups
  (length(unique(snp_group_chr_matrix_wcounts_filtered_nonzero$Res_CHR_Group)) - length(unique(prev$Res_CHR_Group))) + number_filtered_groups + total_number_groups_zero_frac == total_number_groups
  
  
  
  ######## ... TUNING ############
  
  # Tuning works as follows: We know the "fraction" for each group for each gamma. 
  # The "e-value" is just the multiplier "multiplied" by the fraction, so can keep the fraction 
  # and just multiply by all of the different "multipliers"
  
  # first, need unique fractions (not multiple times per group, but only one fraction or e-value per group) 
  unique_filtered_fracs <- filtered_fracs %>% dplyr::select(Res_CHR_Group, fracs) %>% unique()
  
  # Make sure seed is set
  set.seed(2022)
  # generate random matrix consisting of 10000+ rows with random integers summing up to the total_number_groups
  multiplier_matrix <- integer_matrix_fixed_rowsum_UKB(7, total_number_groups, l = 10000, add_extremes = TRUE)
  
  rejections_multiplier_matrix <- c()
  
  ######## ...calculate number of rejections for every row of the multiplier matrix 
  for(m in 1:nrow(multiplier_matrix)) {
    
    if(m %% 500 == 0) {
      print(m)
    }
    
    # extract the row of the multiplier matrix
    multiplier_matrix_row <- data.frame(Resolution = seq(0, 6), Multiplier = multiplier_matrix[m, ])
    
    # for each unique res_chr_group-fraction combination, merger in the corresponding multiplier 
    # so that we can multiply the fraction by the multiplier
    filtered_fracs_merged_mult <- left_join(filtered_fracs, multiplier_matrix_row, by = c("Resolution"))  
    
    # calculate the e-value with the new multiplier
    unique_filtered_fracs <- filtered_fracs_merged_mult %>% 
      mutate(eval_multiplier = fracs*Multiplier) %>%
      dplyr::select(Res_CHR_Group, eval_multiplier) %>% unique()
    
    
    ########### ... RUN e-BH on the newly defined e-values #########
    
    filtered_eval_mult_unique <- unique_filtered_fracs$eval_multiplier
    filtered_fracs_ord <- order(filtered_eval_mult_unique, decreasing = TRUE)
    filtered_eval_mult_unique <- sort(filtered_eval_mult_unique, decreasing = TRUE)
    alpha_mod <- alpha*number_filtered_groups/total_number_groups
    comp <- filtered_eval_mult_unique >= (number_filtered_groups / (alpha_mod * (1:length(filtered_eval_mult_unique)))) 
    
    
    id <- max(which(comp > 0))
    if(id > 0){
      rej <- filtered_fracs_ord[1:id]
    }else{
      rej<- NULL
    }
    
    if(!is.infinite(id)) {
      rejections_multiplier_matrix[m] <- length(rej)
    } else {
      rejections_multiplier_matrix[m] <- 0
    }
    
  }
  
  
  
  ######### ...choose multiplier corresponding to max index ######
  
  # Choose the LAST index (in case multiple  multipliers give the same number of rejections, want the last index 
  # because the last indices in the multiplier matrix correspond to uniform / extreme cases)
  last_max_index <- length(rejections_multiplier_matrix) - which.max(rev(rejections_multiplier_matrix)) + 1
  
  total_rejections_filter_gamma_frac[1, counter]  <- rejections_multiplier_matrix[last_max_index]
  total_rejections_filter_gamma_unif[1, counter] <- rejections_multiplier_matrix[10001] # uniform multiplier index, always

  print("Number of rejections uniform multiplier")
  print(rejections_multiplier_matrix[10001])
    
  print("Number of rejections best multiplier")
  print(rejections_multiplier_matrix[last_max_index])
  
  # save best multiplier 
  best_multiplier_matrix[, counter] <- multiplier_matrix[last_max_index,]
  
  print("Corresponding multiplier")
  print(multiplier_matrix[last_max_index,])
  
  
  ######### ... number of rejections for original e-values ########### 

  # first, need unique e-values (not multiple times per group, but only one per group) 
  unique_filtered_evals <- filtered_fracs %>% dplyr::select(Res_CHR_Group, evals) %>% unique()
  
  filtered_evals_unique <- unique_filtered_evals$evals
  filtered_evals_ord <- order(filtered_evals_unique, decreasing = TRUE)
  filtered_evals_unique <- sort(filtered_evals_unique, decreasing = TRUE)
  alpha_mod <- alpha*number_filtered_groups/total_number_groups
  comp <- filtered_evals_unique >= (number_filtered_groups / (alpha_mod * (1:length(filtered_evals_unique)))) # 7 resolutions in total
  
  
  id <- max(which(comp > 0))
  if(id > 0){
    rej <- filtered_evals_ord[1:id]
  }else{
    rej<- NULL
  }
  

  # which groups were rejected
  rejected_group_evals <- unique_filtered_evals[rej, ]
  rejected_Res_CHR_Group <- rejected_group_evals$Res_CHR_Group
  filtered_rejected <- rejected_group_evals$evals

  rejected_ones_with_info <- filtered_fracs %>% dplyr::filter(Res_CHR_Group %in% rejected_Res_CHR_Group)
  
  # check whether there are any duplicate SNPs
  length(rejected_ones_with_info$SNP) == length(unique(rejected_ones_with_info$SNP))
  
  
  if(!is.infinite(id)) {
    total_rejections_filter_gamma_eval[1, counter] <- length(unique(rej))
  }

  
}

saveRDS(total_rejections_filter_gamma_unif,  paste0(my_out_dir, "RDS/", outcome_string, "_", population, "_fast_total_rejections_filter_gamma_unif_red_f", alpha))
saveRDS(total_rejections_filter_gamma_frac,  paste0(my_out_dir, "RDS/",outcome_string, "_", population, "_fast_total_rejections_filter_gamma_red_f", alpha))
saveRDS(total_rejections_filter_gamma_eval,  paste0(my_out_dir, "RDS/",outcome_string, "_", population, "_fast_evals_total_rejections_filter_gamma_red_f", alpha))   
saveRDS(gamma_matrix,  paste0(my_out_dir, "RDS/",outcome_string, "_", population, "gamma_matrix_red_f", alpha))  
saveRDS(best_multiplier_matrix,  paste0(my_out_dir, "RDS/",outcome_string, "_", population, "best_multiplier_matrix_red_f", alpha))  

