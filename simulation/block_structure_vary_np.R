#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

######################################################################

# Run simulation with block-structured correlation matrix #

######################################################################

############ READ ARGUMENTS & SET BASE PARAMS ############

np_ratio <- as.numeric(args[1])
if(is.na(np_ratio)) np_ratio = 1

sparsity <- as.numeric(args[2])
if(is.na(sparsity)) sparsity = 0.1

nrep <- as.numeric(args[3])
if(is.na(nrep)) nrep <- 5

rho <- as.numeric(args[4])
if(is.na(rho)) rho <- 0.8

M <- as.numeric(args[5])
if(is.na(M)) M = 2

p <- as.numeric(args[6])
if(is.na(p)) p = 1000 
#2500

block_sigma_ind <- as.numeric(args[7])
if(is.na(block_sigma_ind)) block_sigma_ind = 1

if(block_sigma_ind == 1) {
  block_sigma = TRUE
} else {
  block_sigma = FALSE
}

gamma_tune_threshold <- as.numeric(args[8])
if(is.na(gamma_tune_threshold)) gamma_tune_threshold = 0.5

sd_dist_beta <- as.numeric(args[9])
if(is.na(sd_dist_beta)) sd_dist_beta = 0.2


block_sigma_type_indicator <- as.numeric(args[10])
if(is.na(block_sigma_type_indicator)) block_sigma_type_indicator = 1


# simulation setup - set paths
if(dir.exists("")){
  mydir <- ""
  setwd(mydir)
  
  my_out_dir <- "" 
  my_julia_path <- ""
} 

## The directory to save the results
save_dir <- sprintf(paste0(mydir, "output/simulation"))
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}



print(paste0("np_ratio is", np_ratio))
print(paste0("nrep is", nrep))
print(paste0("rho is", rho))
print(paste0("M is", M))
print(paste0("p is", p))
print(paste0("block_sigma_ind is", block_sigma_ind))
print(paste0("sd_dist_beta is", sd_dist_beta))


# LOAD ALL FUNCTIONS 
myseed = 2022

set.seed(myseed)

source(paste0(mydir, "utils/load_aux.R"))

# tell JuliaCall where is the Julia executable (the path can be found by typing `which julia` on terminal)
julia <- julia_setup(JULIA_HOME = my_julia_path)

# load Knockoffs.jl package
julia_library("Knockoffs")

# set parameters of problem 
n = round(np_ratio * p)

print("N is")
print(n)

# number of outcomes
Q = 1
# percentage of nonzero snps that are the same across outcomes
overlap_pct = 1
# X follows AR(k) process
# set specific parameters for simulating X

# knockoff parameters
mu = rep(0, p)
partial_u = 1
method = "maxent"

# FDR control parameters
alpha = 0.2


############### FIX GROUP STRUCTURE ############


if(M == 1) {
  group_size_with_individual <- c(1)
  
} 

if(M == 2) {
  group_size_with_individual <- c(1, 5)
} 

if(M == 3) {
  group_size_with_individual <- c(1, 5, 10)
}

if(M == 4) {
  group_size_with_individual <- c(1, 5, 10, 20)
}

if(block_sigma_type_indicator == 1) {
  group_size_for_block_sigma = group_size_with_individual[2]
  block_size_vector = rep(group_size_for_block_sigma, p / group_size_for_block_sigma) 
}

if(block_sigma_type_indicator == 2) {
  
  vector_group_10 <- rep(10, floor(0.5*p/10))
  vector_group_5 <- rep(5, floor(0.5*p/5))
  
  block_size_vector <- sample(c(vector_group_10, vector_group_5))
  
}




number_cutoffs = length(group_size_with_individual)
method <- "maxent"


group_struct_res <- generate_group_structure(group_size_with_individual, p, M) 

# groups is a data frame of dimension p x M 
# first column has group number in individual layer 
# second column has group number in second layer 
# ... and so on
groups_combined <- group_struct_res[[1]]

if(M == 2) {
  colnames(groups_combined) <- c("individual", "layer2") 
} else if(M == 3) {
  colnames(groups_combined) <- c("individual", "layer2", "layer3") 
} else {
  colnames(groups_combined) <- c("individual", "layer2", "layer3", "layer4")
}





# group_candidates_list is a list whith length == total number of groups (i.e. all individuals + all groups in other layers) 
# each list element contains the number of the individual elements that are part of each group
cand_groups <- group_struct_res[[2]]

group_data_frame <- group_struct_res[[3]]



number_groups_in_res <- c() 
for(c in 1:M) {
  number_groups_in_res[c] <- p / group_size_with_individual[c]
}
############### FIX BETA ######################


# CREATE BETA
beta_res <- random_beta_generation(p = p, Q = Q, sparsity = sparsity, overlap_pct = overlap_pct, myseed = myseed)
beta <- beta_res$beta_ind_random_sign
beta_ind_nonzero <- beta_res$nonzero_matrix

nonzero_per_q <- beta_res$nonzero_per_q
nonzero_any_indiv <- beta_res$nonzero_any_indiv
print("NUMBER NONZERO")
print(length(nonzero_any_indiv))

number_cutoffs = M

nonzero_group_per_layer = sapply(1:M, function(m)(unique(groups_combined[nonzero_any_indiv,m]))) 


############### SET GAMMA ###############

# this is used in the paper: gamma = alpha / 2 = 0.2 / 2 = 0.1
if(gamma_tune_threshold == 1) {
  gammas <- rep(alpha/2, M)
} else {
  print("Need to specify gamma.")
}

############### BEGIN REP SIM LOOP #############


all_res_by_resolutions <- c()
all_res_across_resolutions <- c()
my_snr_df <- c()

for(r in 1:nrep) {
  
  seedB <- myseed + r
  set.seed(seedB)
  
  cat(sprintf("Running the %d-th rep.\n", r)) 
  
  # define objects to hold simulation results
  groups_list <- c()
  W_list <- c()
  fracs_list <- list()
  rejected_knockoff <- list()
  nonzero_group_per_layer <- list()
  nonzero_res_group_per_layer <-list()
  
  knockoff_fdr_per_layer <- matrix(NA, nrow = 1, ncol = number_cutoffs)
  knockoff_power_per_layer <- matrix(NA, nrow = 1, ncol = number_cutoffs)
  knockoff_true_power_per_layer <- matrix(NA, nrow = 1, ncol = number_cutoffs)
  knockoff_snps_implicated_per_layer <- matrix(NA, nrow = 1, ncol = number_cutoffs)
  knockoff_frac_indiv_correct_per_layer <- matrix(NA, nrow = 1, ncol = number_cutoffs)
  
  kelp_fdr = matrix(NA, nrow = 1, ncol = 1)
  kelp_power = matrix(NA, nrow = 1, ncol = 1)
  kelp_snps_implicated = matrix(NA, nrow = 1, ncol = 1)
  kelp_frac_correct = matrix(NA, nrow = 1, ncol = 1)
  
  data_res <- generate_data(Sigma_Dirichlet = FALSE, 
                            block_sigma = block_sigma, # USE BLOCK STRUCTURE
                            distribution_beta = TRUE, # USE DISTRIBUTION BETA
                            binary_y = FALSE,
                            beta = beta,
                            n = n,
                            p = p,
                            rho = rho,
                            Q = Q, 
                            mean_dist_beta = 0,
                            sd_dist_beta = sd_dist_beta,
                            myseed = seedB, 
                            groups_combined, 
                            group_size_with_individual = group_size_with_individual, 
                            block_size_vector = block_size_vector)
  
  
  Y <- data_res$Y
  X <- data_res$X
  
  etas <- data_res$etas
  Sigma <- data_res$Sigma
  
  my_snr <- data_res$my_snr
  my_snr_alt <- data_res$my_snr_alt
  
  cand_groups <- list()
  
  counter_layer = 1
  
  
  W_list <- list()
  knockoff_fraction_per_layer <- list()
  knockoff_eval_per_layer <- list()
  selected_knockoff_per_layer <- list()
  
  
  ############### ... BEGIN CUTOFF LOOP #############
  for(c in 1:M) {
    
    
    selected_knockoff_groups_per_q <- list()
    knockoff_fraction_per_layer[[c]] <- list()
    W_list[[c]] <- list()
    
    cat(paste0("Processing cutoff ", c))
    
    
    groups_df = data.frame(groups_combined[, c(1, c)])
    colnames(groups_df) <- c("individual", "groups")
    number_unique_groups <- length(unique(groups_df$groups))
    groups = groups_df$groups
    
    nonzero_group_numbers_per_layer <- groups_df %>% dplyr::filter(individual %in% nonzero_any_indiv) %>% pull(groups) %>% unique()
    nonzero_group_per_layer[[counter_layer]] <- nonzero_group_numbers_per_layer
    nonzero_res_group_per_layer[[counter_layer]] <- paste0("res_", c, "_group_",nonzero_group_numbers_per_layer)
    
    
    # create cand group list 
    temp <- split(groups_df[, 1],groups_df[, 2])
    group_number <- paste0("group_", names(temp))
    
    cand_groups_layer <- mapply(list, SNP = temp, group = group_number, Resolution = c, SIMPLIFY=FALSE)
    
    cand_groups_layer <- lapply(cand_groups_layer,
                                function(x) c(x, Res_Group =paste0("res_", x$Resolution, "_", x$group)))
    
    cand_groups_layer <- unname(cand_groups_layer)
    
    cand_groups <- c(cand_groups, cand_groups_layer)
    
    fracs_per_q <- matrix(NA, nrow = length(unique(groups)), ncol = Q)
    
    
    
    ############### ...... BEGIN OUTCOME LOOP #############
    
    for(q in 1:Q) {
      
      y = as.matrix(Y[, q])
      
      if((c == 1) | (number_unique_groups == p)) {
        result <- run_ko_lasso_with_julia(y = y, X = X, Sigma = Sigma, mu = mu, method = method, m = 1, alpha = alpha)
        
      } else {
        # call Julia solver
        result <- run_group_ko_lasso_with_julia(y = y, X = X, Sigma = Sigma, mu = mu,
                                                method = method, m = 1, groups = groups, alpha = alpha)
      }
      
      
      
      # get rejected groups by the knockoff 
      selected_knockoff_groups_per_q[[q]] <- result$selected[1]
      
      print(length(result$selected[1]))
      
      # get fractions for e-values
      W_list[[c]][[q]] <- result$W
      
      
      tau_knockoff <- knockoff.threshold(W_list[[c]][[q]], fdr=alpha, offset=1)
      knockoff_fraction_per_layer[[c]][[q]] <- (W_list[[c]][[q]] >= tau_knockoff) / (1 + sum(W_list[[c]][[q]] <= -tau_knockoff))
      knockoff_eval_per_layer <- (p / group_size_with_individual[c]) * (W_list[[c]][[q]] >= tau_knockoff) / (1 + sum(W_list[[c]][[q]] <= -tau_knockoff))
      
    }
    
    selected_knockoff_per_layer[[c]] <- selected_knockoff_groups_per_q
    
    
    ############### ...... END OUTCOME LOOP + COMBINE RES #############
    
    
    ############## ................ knockoff per layer #########################
    
    
    # reject the knockoff if it is rejected for any outcome
    rejected_knockoff <- unique(unlist(selected_knockoff_groups_per_q))
    
    if(length(rejected_knockoff) > 0) {
      knockoff_fdr_per_layer[1, counter_layer] = length(setdiff(rejected_knockoff, nonzero_group_per_layer[[counter_layer]])) / length(rejected_knockoff)
      
      knockoff_power_per_layer[1, counter_layer] = (length(intersect(rejected_knockoff,  nonzero_group_per_layer[[counter_layer]])) / group_size_with_individual[counter_layer]) / length(nonzero_any_indiv)
      knockoff_true_power_per_layer[1, counter_layer] = length(intersect(rejected_knockoff,  nonzero_group_per_layer[[counter_layer]])) / length(nonzero_group_per_layer[[counter_layer]])
      
      
      
      # calculate number of SNPs implicated and fraction correct rejections
      knockoff_snps_implicated_per_layer[1, counter_layer] = nrow(groups_df %>% dplyr::filter(groups %in% rejected_knockoff))
      number_correct_rej_individual <- nrow(groups_df %>% dplyr::filter(groups %in% rejected_knockoff) %>% dplyr::filter(individual %in% nonzero_any_indiv))
      knockoff_frac_indiv_correct_per_layer[1, counter_layer] = number_correct_rej_individual / length(nonzero_any_indiv)
      
      print("SNPs implicated, Frac correct")
      print(knockoff_snps_implicated_per_layer)
      print(knockoff_frac_indiv_correct_per_layer)
    } else {
      knockoff_fdr_per_layer[1, counter_layer]  = 0
      knockoff_power_per_layer[1, counter_layer] = 0
      knockoff_true_power_per_layer[1, counter_layer] = 0
      knockoff_frac_indiv_correct_per_layer[1, counter_layer] = 0
      knockoff_snps_implicated_per_layer[1, counter_layer] = 0
    }
    
    
    
    
    # increase counter for next layer
    counter_layer = counter_layer + 1
  }
  
  
  ############## ................ knockoff outer  #########################
  
  knockoff_rejections_info <- c() 
  
  for(c in 1:M) {
    knockoff_rejections_info_m <- group_data_frame %>% filter(Resolution == c, 
                                                              group_number %in% selected_knockoff_per_layer[[c]][[1]])
    knockoff_rejections_info <- rbind(knockoff_rejections_info, knockoff_rejections_info_m)
  }
  
  
  outer_nodes_knockoff <- filter_outer_node(knockoff_rejections_info,
                                            resolutions = seq(1, M),
                                            group_identifiyer = "res_Group")
  
  selected_outer_knockoff_groups <- outer_nodes_knockoff$res_Group 
  
  # calculate power and fdr for knockoff outer 
  if(length(selected_outer_knockoff_groups) > 0) {
    
    # get power, fdr and frac correct rejections
    knockoff_outer_fdr = length(setdiff(selected_outer_knockoff_groups, unlist(nonzero_res_group_per_layer))) / length(selected_outer_knockoff_groups)
    
    
    true_group_rejections <- intersect(selected_outer_knockoff_groups, unlist(nonzero_res_group_per_layer))
    # get the size for each of these true rejections
    group_data_frame %>% 
      dplyr::filter(res_Group %in% unlist(nonzero_res_group_per_layer)) %>%
      dplyr::filter(res_Group %in% true_group_rejections) %>% 
      group_by(res_Group) %>% 
      count() %>% 
      mutate(inverse_size = 1 / n) -> inverse_size_true_knockoff_outer_rejections
    
    group_data_frame %>% 
      dplyr::filter(res_Group %in% selected_outer_knockoff_groups) %>% 
      pull(SNP) %>% 
      unique() -> knockoff_outer_rejected_snps
    
    # sum of inverse size of ture rejections divided by total number of nonzero snps
    knockoff_outer_power <- sum(inverse_size_true_knockoff_outer_rejections$inverse_size) / length(nonzero_group_per_layer[[1]])
    
    knockoff_outer_snps_implicated <- length(knockoff_outer_rejected_snps)
    
    knockoff_outer_frac_correct = length(intersect(knockoff_outer_rejected_snps, nonzero_group_per_layer[[1]])) / length(nonzero_group_per_layer[[1]])
    
  } else {
    knockoff_outer_fdr = 0
    knockoff_outer_power = 0
    knockoff_outer_frac_correct = 0
    knockoff_outer_snps_implicated = 0
    
    ebh_knockoff_outer_fdr = 0
    ebh_knockoff_outer_power = 0
    ebh_knockoff_outer_snps_implicated = 0
    ebh_knockoff_outer_frac_correct = 0
  }
  
  
  ####################### .........................  + e-BH ########################
  
  
  if(length(selected_outer_knockoff_groups) > 0 & !(is.na(selected_outer_knockoff_groups)[1])) {
    
    # create structure to hold outer nodes with fractions
    outer_nodes_knockoff_w_fracs <- group_data_frame %>% 
      dplyr::filter(res_Group %in% outer_nodes_knockoff$res_Group) %>%
      dplyr::select(group_number, Resolution, res_Group) %>%
      unique() %>%
      mutate(fracs_knockoff = NA)
    
    # add knockoff fractions 
    for(c in 1:M) {
      group_numbers <- outer_nodes_knockoff_w_fracs %>% filter(Resolution == c) %>% pull(group_number) %>% unique()
      
      for(g in group_numbers) {
        outer_nodes_knockoff_w_fracs[outer_nodes_knockoff_w_fracs$Resolution == c & outer_nodes_knockoff_w_fracs$group_number == g, ]$fracs_knockoff <- knockoff_fraction_per_layer[[c]][[1]][g]
      }
    }
    
    
    total_number_groups <- length(cand_groups)
    number_groups_outer <- length(unique(outer_nodes_knockoff_w_fracs$res_Group))
    alpha_mod <- alpha * number_groups_outer / total_number_groups
    # e-values based on the knockoff are defined as the original e-values!
    outer_nodes_knockoff_w_fracs <- outer_nodes_knockoff_w_fracs %>% 
      mutate(eval_multiplier = number_groups_in_res[Resolution]) %>% 
      mutate(eval_calc = eval_multiplier * fracs_knockoff)
    
    ebh_outer_nodes_rej_indices <- run_eBH(outer_nodes_knockoff_w_fracs$eval_calc, alpha_mod, 
                                           number_groups = number_groups_outer, use_fractions = FALSE )
    
    
    # calculate power, fdr and snps implicated
    
    if(!is.na(ebh_outer_nodes_rej_indices[1])) {
      
      ebh_knockoff_outer_rejected_res_group <- outer_nodes_knockoff_w_fracs[ebh_outer_nodes_rej_indices,]$res_Group
      
      # get power, fdr and frac correct rejections
      ebh_knockoff_outer_fdr = length(setdiff(ebh_knockoff_outer_rejected_res_group, unlist(nonzero_res_group_per_layer))) / length(ebh_knockoff_outer_rejected_res_group)
      
      
      true_group_rejections_knockoff_ebh <- intersect(ebh_knockoff_outer_rejected_res_group, unlist(nonzero_res_group_per_layer))
      # get the size for each of these true rejections
      group_data_frame %>% 
        dplyr::filter(res_Group %in% unlist(nonzero_res_group_per_layer)) %>%
        dplyr::filter(res_Group %in% true_group_rejections_knockoff_ebh) %>% 
        group_by(res_Group) %>% 
        count() %>% 
        mutate(inverse_size = 1 / n) -> inverse_size_true_ebh_knockoff_outer_rejections
      
      group_data_frame %>% 
        dplyr::filter(res_Group %in% ebh_knockoff_outer_rejected_res_group) %>% 
        pull(SNP) %>% 
        unique() -> ebh_knockoff_outer_rejected_snps
      
      # sum of inverse size of ture rejections divided by total number of nonzero snps
      ebh_knockoff_outer_power <- sum(inverse_size_true_ebh_knockoff_outer_rejections$inverse_size) / length(nonzero_group_per_layer[[1]])
      
      ebh_knockoff_outer_snps_implicated <- length(ebh_knockoff_outer_rejected_snps)
      
      ebh_knockoff_outer_frac_correct = length(intersect(ebh_knockoff_outer_rejected_snps, nonzero_group_per_layer[[1]])) / length(nonzero_group_per_layer[[1]])
      
      
    } else {
      ebh_knockoff_outer_fdr = 0
      ebh_knockoff_outer_power = 0
      ebh_knockoff_outer_snps_implicated = 0
      ebh_knockoff_outer_frac_correct = 0
    }
    
  }
  
  
  
  ############## ................ kelp #########################
  # begin loop to construct e-values
  
  for(c in 1:M) {
    
    fracs_per_q <- matrix(NA, nrow = length(unique(groups_combined[, c(c)])), ncol = Q)
    
    for(q in 1:Q) {
      
      tau_gamma <- knockoff.threshold(W_list[[c]][[q]], fdr=gammas[c], offset=1)
      print(tau_gamma)
      fracs_per_q[, q] <- (W_list[[c]][[q]] >= tau_gamma) / (1 + sum(W_list[[c]][[q]] <= -tau_gamma))
      
      
    }
    
    fracs_list[[c]] <- rowMeans(fracs_per_q)
    
  }
  
  
  # add e-values to cand list
  cand_groups_kelp <- cand_groups
  
  for(g in 1:length(cand_groups)) {
    cand_groups_kelp[[g]]$fracs <- unlist(fracs_list)[g]
  }
  
  # get kelp rejections    
  # filter out groups with zero e-values (or fractions)
  cand_groups_kelp_filtered <- cand_groups_kelp[sapply(cand_groups_kelp, function(x) x$fracs > 0)]
  
  # get individual "snps" in the non-zero groups
  individual_snps_candidates <-  unique(unlist(sapply(cand_groups_kelp_filtered, `[[`, "SNP")))
  
  # number groups in optimization procedure
  n_groups_filtered <- length(cand_groups_kelp_filtered)
  
  # number of unique snps 
  n_individual_hypotheses_filtered <- length(individual_snps_candidates)
  
  # run optimization problem
  if(n_groups_filtered > 0) {
    cat("Running kelp ...")
    res <- kelp(group_candidates_list = cand_groups_kelp_filtered, 
                individual_snps_candidates = individual_snps_candidates, 
                single_u = TRUE,
                alpha = alpha, 
                verbose = F,
                weighted = T, 
                fractions_only = TRUE,
                M = number_cutoffs, 
                n_individual_hypotheses_filtered = n_individual_hypotheses_filtered, 
                n_groups_filtered = n_groups_filtered, 
                total_groups = total_number_groups, 
                Q = Q,
                u_seq = partial_u) # this automatically adjust alpha to the required level
    cat("done\n.")
    kelp_detections <- res[[1]]
    
    # get selected groups
    selected_kelp_groups <- sapply(kelp_detections, `[[`, "Res_Group")
    
    selected_kelp_fracs <- sapply(kelp_detections, `[[`, "fracs")
    included_kelp_SNPs <- unlist(sapply(kelp_detections, `[[`, "SNP"))
    
  } else {
    selected_kelp_groups <- NULL
  }
  
  # get power and fdr for kelp  # check if self consistent
  
  if(length(selected_kelp_groups) > 0) {
    # get power, fdr and frac correct rejections
    kelp_fdr[1, 1] = length(setdiff(selected_kelp_groups, unlist(nonzero_res_group_per_layer))) / length(selected_kelp_groups)
    
    
    true_group_rejections <- intersect(selected_kelp_groups, unlist(nonzero_res_group_per_layer))
    kelp_rejections_true_filtered <- kelp_detections[sapply(kelp_detections, function(x) x$Res_Group %in% true_group_rejections)]
    
    
    kelp_power[1, 1] <- sum(sapply(kelp_rejections_true_filtered, inverse_size)) / length(nonzero_any_indiv)
    
    kelp_snps_implicated[1, 1] <- sum(sapply(kelp_detections, function(x) length(x$SNP)))
    kelp_frac_correct[1, 1] = length(intersect(included_kelp_SNPs, nonzero_any_indiv)) / length(nonzero_any_indiv)
  } else {
    kelp_fdr[1, 1]= 0
    kelp_power[1, 1] = 0
    kelp_snps_implicated[1, 1] = 0
    kelp_frac_correct[1, 1] = 0
  }
  
  
  print(kelp_power)
  print(kelp_snps_implicated)
  print(kelp_frac_correct)
  
  
  ########### ... END CUTOFF LOOP #################
  summary_by_resolutions = data.frame(level = rep("by_resolution", number_cutoffs),
                                      method = c(rep("knockoff", number_cutoffs)),
                                      Resolution = seq(1:M),
                                      fdp = c(knockoff_fdr_per_layer),
                                      power = c(knockoff_power_per_layer),
                                      true_power = c(knockoff_true_power_per_layer),
                                      snps_implicated = c(knockoff_snps_implicated_per_layer), 
                                      frac_correct_rejections = c(knockoff_frac_indiv_correct_per_layer),
                                      np_ratio = rep(np_ratio, number_cutoffs), 
                                      seed = rep(seedB, number_cutoffs), 
                                      Q = rep(Q, number_cutoffs), 
                                      overlap_pct = rep(overlap_pct, number_cutoffs), 
                                      alpha = rep(alpha, number_cutoffs))
  
  
  summary_across_resolutions = data.frame(level = rep("across_resolution", 3),
                                          method = c("kelp", "knockoff_outer", "knockoff_outer_ebh"),
                                          fdp = c(kelp_fdr, knockoff_outer_fdr, ebh_knockoff_outer_fdr),
                                          power = c(kelp_power, knockoff_outer_power, ebh_knockoff_outer_power),
                                          snps_implicated = c(kelp_snps_implicated, knockoff_outer_snps_implicated, ebh_knockoff_outer_snps_implicated), 
                                          frac_correct_rejections = c(kelp_frac_correct, knockoff_outer_frac_correct, ebh_knockoff_outer_frac_correct),
                                          np_ratio = rep(np_ratio, 3), 
                                          seed = rep(seedB, 3), 
                                          Q = rep(Q, 3), 
                                          overlap_pct = rep(overlap_pct, 3), 
                                          alpha = rep(alpha, 3), 
                                          gamma = rep(paste0(gammas, collapse = ""), 3))
  
  
  summary_my_snr = data.frame(my_snr = my_snr, 
                              my_snr_alt = my_snr_alt)
  
  
  all_res_by_resolutions <- rbind(all_res_by_resolutions, summary_by_resolutions)
  all_res_across_resolutions <- rbind(all_res_across_resolutions, summary_across_resolutions)
  my_snr_df <- rbind(my_snr_df, summary_my_snr)
  
}


############### SAVE RESULTS #############

print(all_res_by_resolutions %>% group_by(method, Resolution, np_ratio) %>% summarize(mean_fdr = mean(fdp, na.rm = TRUE), 
                                                                                 mean_power = mean(power, na.rm = TRUE), 
                                                                                 mean_true_power = mean(true_power, na.rm = TRUE), 
                                                                                 mean_snps_implicated = mean(snps_implicated), 
                                                                                 mean_frac_correct = mean(frac_correct_rejections)))

print(all_res_across_resolutions %>% group_by(method, np_ratio) %>% summarize(mean_fdr = mean(fdp, na.rm = TRUE), 
                                                                         mean_power = mean(power, na.rm = TRUE), 
                                                                         mean_snps_implicated = mean(snps_implicated), 
                                                                         mean_frac_correct = mean(frac_correct_rejections)))




## Save the outcomes
out_dir <- sprintf("%s/np_block_beta_dist_np_ratio_%dM_g_res_by_Q_%d_np_ratio_%.1f_nrep_%d_ol_%.1f_g%s_al%.4f_u%d_ncuts_%d_%s_r%.3f_p%d_s%.3f_gtt%.3f_sdb%.3f.csv", save_dir, 
                   M, Q,
                   np_ratio, nrep, overlap_pct, paste0(gammas, collapse = "_"), 
                   alpha, partial_u, number_cutoffs, paste0(group_size_with_individual, collapse = "_"), 
                   rho, p, sparsity, gamma_tune_threshold, sd_dist_beta)
write_csv(all_res_by_resolutions, out_dir)

out_dir <- sprintf("%s/np_block_beta_dist_res_across_block_np_ratio_%dM_g_res_by_Q_%d_np_ratio_%.1f_nrep_%d_ol_%.1f_g%s_al%.4f_u%d_ncuts_%d_%s_r%.3f_p%d_s%.3f_gtt%.3f_sdb%.3f.csv", save_dir, 
                   M, Q,
                   np_ratio, nrep, overlap_pct, paste0(gammas, collapse = "_"), 
                   alpha, partial_u, number_cutoffs, paste0(group_size_with_individual, collapse = "_"), 
                   rho, p, sparsity, gamma_tune_threshold, sd_dist_beta)
write_csv(all_res_across_resolutions, out_dir)

out_dir <- sprintf("%s/np_block_beta_dist_snr_np_ratio_%dM_g_res_by_Q_%d_np_ratio_%.1f_nrep_%d_ol_%.1f_g%s_al%.4f_u%d_ncuts_%d_%s_r%.3f_p%d_s%.3f_gtt%.3f_sdb%.3f.csv", save_dir, 
                   M, Q,
                   np_ratio, nrep, overlap_pct, paste0(gammas, collapse = "_"), 
                   alpha, partial_u, number_cutoffs, paste0(group_size_with_individual, collapse = "_"), 
                   rho, p, sparsity, gamma_tune_threshold, sd_dist_beta)
write_csv(data.frame(my_snr_df), out_dir)




