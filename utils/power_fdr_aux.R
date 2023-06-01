



################### HELPER FUNCTIONS FOR POWER & FDP  ###############################


# function for power and fdr calculation given list with rejections

# INPUT: 

# rejections_list: list containing rejections at each level of Resolution 
# nonzero_groups_list: list containing the true non-zero groups at each level of Resolution 
# M: number of resolutions (integer)
# by_resolution: indicator for whether power / fdr should be calculated separately by resolution or overall
# group_size_with_individual: vector containing group sizes, including the individual layer (e.g. c(1, 5, 10)) for group sizes 1 (individual), 5 and 10
# alpha: target FDR control
# group_data_frame: data frame containing group information

# OUTPUT: List containing calculated power, fdr, snps implicated and fraction correct


# DOESN"T SEEM LIKE IM USING THIS FUNCTION 

power_fdp_calc <- function(rejections_list, nonzero_groups_list, M, by_resolution, 
                           group_size_with_individual, alpha, group_data_frame) {
  
  
  
  
  if(length(rejections_list) > 0) {
    
    if(by_resolution) {
      
      fdp <- matrix(NA, nrow = 1, ncol = M)
      power <- matrix(NA, nrow = 1, ncol = M)
      snps_implicated <- matrix(NA, nrow = 1, ncol = M)
      frac_correct_rejections <- matrix(NA, nrow = 1, ncol = M)
      
      
      for(m in 1:M) {
        
        
        if(!((length(rejections_list[[m]])==1 && is.na(rejections_list[[m]])))) {
          fdp[, m] <-  length(setdiff(rejections_list[[m]],
                                      nonzero_groups_list[[m]]))/length(rejections_list[[m]]) 
          
          power[, m] <-  length(intersect(rejections_list[[m]],
                                          nonzero_groups_list[[m]]))/length(nonzero_groups_list[[m]]) 
          
          # calculate snps implicitly rejected and pct of individual signals correctly rejected
          
          group_data_frame_rejected <- group_data_frame %>% dplyr::filter(Resolution == m, SNP %in% rejections_list[[m]])
          
          snps_implicated[, m] <- nrow(group_data_frame_rejected)
          frac_correct_rejections[, m] <- length(intersect(group_data_frame_rejected$SNP, nonzero_groups_list[[1]])) / length(nonzero_groups_list[[1]])
        }
        
        
      }
      
    } else {
      
      fdp = sum(sapply(1:M, function(m)(length(setdiff(rejections_list[[m]], nonzero_groups_list[[m]]))))) / sum(sapply(rejections_list, length))
      
      # power = weight-adjusted number of rejections
      number_rejections_by_res <- sapply(rejections_list, length)
      weight = 1 / group_size_with_individual
      power = sum(number_rejections_by_res * weight)
      power = power / length(nonzero_groups_list[[1]])
      
      # calculate snps implicitly rejected and pct of individual signals correctly rejected
      group_data_frame_rejected <- c()
      for(m in 1:M) {
        
        group_data_frame_rejected <- rbind(group_data_frame_rejected, 
                                           group_data_frame %>% 
                                             dplyr::filter(Resolution == m, SNP %in% rejections_list[[m]]))
        
        
      }
      
      
      snps_implicated <- nrow(group_data_frame_rejected)
      frac_correct_rejections <- length(intersect(group_data_frame_rejected$SNP, nonzero_groups_list[[1]])) / length(nonzero_groups_list[[1]])
      
      
      
    }
    
  }
  
  
  fdp[is.na(fdp)] <- 0
  power[is.na(power)] <- 0
  
  return(list(fdp = fdp, power = power, 
              snps_implicated = snps_implicated,
              frac_correct_rejections = frac_correct_rejections))
  
}


# function for power and fdr calculation given list with rejections for hierarchically structured outcomes

# INPUT: 
# knockoff_rejected: data frame containig the knockoff rejections (one row = one rejection)
# nonzero_res_groups: vector containing all nonzero groups 
# nonzero_groups_per_q: List containing nonzero groups, each sub-list belongs to one outcome

# OUTPUT: list containing power, fdp, fraction correct rejected, number of implicated outcomes

power_fdp_calc_knockoff_per_resolution_hierarchical <- function(knockoff_rejected, 
                                                                nonzero_snps_per_q_character,
                                                                nonzero_snps_per_q, 
                                                                nonzero_individual_snps) {
  
  # power for knockoff rejections per layer
  # only count correct rejections!
  knockoff_rejected %>% 
    mutate(res_group = paste0("outcome_", outcome, "_SNP_", SNP)) %>% 
    dplyr::filter(res_group %in% nonzero_snps_per_q_character) %>%
    group_by(resolution) %>% 
    summarise(power = sum(level) / length(nonzero_individual_snps)) -> knockoff_power_by_resolution_df
  
  knockoff_rejected %>% 
    mutate(res_group = paste0("outcome_", outcome, "_SNP_", SNP)) %>% 
    dplyr::filter(res_group %in% nonzero_snps_per_q_character) %>%
    group_by(resolution) %>% 
    summarise(frac_correct = n() / length(nonzero_individual_snps)) -> knockoff_frac_correct_by_resolution_df
  
  knockoff_n_outcomes_implicated_df <- knockoff_rejected %>% 
    group_by(resolution) %>% 
    summarise(n_outcomes_implicated = sum(n_outcomes))
  
  
  # force output to have 3 entries
  available_res <- data.frame(resolution = seq(1, 3), power_all_na = rep(NA, 3))
  available_res <- left_join(available_res, knockoff_power_by_resolution_df, by = "resolution")
  available_res <- left_join(available_res, knockoff_frac_correct_by_resolution_df, by = "resolution")
  available_res <- left_join(available_res, knockoff_n_outcomes_implicated_df, by = "resolution")
  
  knockoff_power_by_resolution = available_res$power
  knockoff_frac_correct_by_resolution = available_res$frac_correct
  knockoff_n_outcomes_implicated = available_res$n_outcomes_implicated
  
  # replace NA with 0
  knockoff_frac_correct_by_resolution[is.na(knockoff_frac_correct_by_resolution)] <- 0
  knockoff_power_by_resolution[is.na(knockoff_power_by_resolution)] <- 0
  
  # get power and fdp by level for the knockoff
  knockoff_rejected_by_outcome <- list()
  
  for(q in seq(1, 7)) {
    knockoff_rejected_by_outcome[[q]] <- knockoff_rejected %>% dplyr::filter(outcome == q) %>% pull(SNP)
  }
  
  if(length(knockoff_rejected_by_outcome) > 0) {
    
    fdp_numerators = sapply( seq(1, 7), function(q)(length(setdiff(knockoff_rejected_by_outcome[[q]], nonzero_snps_per_q[[q]]))))
    fdp_denominators = sapply(knockoff_rejected_by_outcome, length)
    
    # resolution 1
    knockoff_fdp_res1 = sum(fdp_numerators[1:4]) / sum(fdp_denominators[1:4])
    knockoff_fdp_res2 = sum(fdp_numerators[5:6]) / sum(fdp_denominators[5:6])
    knockoff_fdp_res3 = sum(fdp_numerators[7]) / sum(fdp_denominators[7])
    
    knockoff_fdp_by_resolution <- c(knockoff_fdp_res1, knockoff_fdp_res2, knockoff_fdp_res3)
    
  } else {
    knockoff_fdp_by_resolution = rep(0, 3)
  }  
  
  return(list(knockoff_fdp_by_resolution = knockoff_fdp_by_resolution, 
              knockoff_power_by_resolution = knockoff_power_by_resolution, 
              knockoff_frac_correct_by_resolution = knockoff_frac_correct_by_resolution, 
              knockoff_n_outcomes_implicated = knockoff_n_outcomes_implicated))
}



# function for calculating power in hierarchical setting (spanning multiple levels of resolution)
# INPUT: 
# rejections: data frame containing the rejections (one row = one rejection)
# nonzero_snps_per_q_character: vector containing all nonzero groups 
# all_outcomes: vector containing the outcome identifications
# nonzero_snps_per_q: List containing nonzero groups, each sub-list belongs to one outcome

# OUTPUT: list containing power, fdp, fraction correct rejected, number of implicated outcomes


power_fdp_calc_hierarchical <- function(nonzero_snps_per_q, all_outcomes, rejections, nonzero_snps_per_q_character, nonzero_individual_snps) {
  
  
  rejections <- rejections %>% mutate(res_group = paste0("outcome_", outcome, "_SNP_", SNP))
  
  n_outcomes_implicated <- sum(rejections$n_outcomes)
  
  true_rejections <- intersect(rejections$res_group, nonzero_snps_per_q_character)
  
  
  rejections_list <- list() 
  
  for(q in all_outcomes) {
    rejections_list[[q]] <- rejections %>% dplyr::filter(outcome == q) %>% pull(SNP)
  }
  
  if(length(rejections_list) > 0) {
    
    fdp = sum(sapply(all_outcomes, function(q)(length(setdiff(rejections_list[[q]], nonzero_snps_per_q[[q]]))))) / sum(sapply(rejections_list, length))
    
    # sum of levels ONLY for true rejections
    power = sum(rejections %>% dplyr::filter(res_group %in% true_rejections) %>% pull(level)) / length(nonzero_individual_snps)
    
    frac_correct_rejections = nrow(rejections %>% dplyr::filter(res_group %in% true_rejections)) / length(nonzero_individual_snps)
    
  } else {
    fdp = 0
    power = 0
    frac_correct_rejections = 0
  }
  
  return(list(fdp = fdp, power = power, frac_correct_rejections = frac_correct_rejections, n_outcomes_implicated = n_outcomes_implicated))
  
  
}

