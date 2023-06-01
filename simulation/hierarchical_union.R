#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)


######################################################################

# Run simulation with hierarchically structured phenotypes #

######################################################################

np_ratio <- as.numeric(args[1])
if(is.na(np_ratio)) np_ratio <- 1

sd_dist_beta <- as.numeric(args[2])
if(is.na(sd_dist_beta)) sd_dist_beta <- 0.5

mean_dist_beta <- as.numeric(args[3])
if(is.na(mean_dist_beta)) mean_dist_beta <- 1

alpha <- as.numeric(args[4])
if(is.na(alpha)) alpha <- 0.2

sparsity <- as.numeric(args[5])
if(is.na(sparsity)) sparsity <- 0.2

overlap_pct <- as.numeric(args[6])
if(is.na(overlap_pct)) overlap_pct <- 0.25

nrep <- as.numeric(args[7])
if(is.na(nrep)) nrep <- 100

p <- as.numeric(args[8])
if(is.na(p)) p <- 1000

proportion_of_ones_desired <- as.numeric(args[9])
if(is.na(proportion_of_ones_desired)) proportion_of_ones_desired <- 0.15


notunegamma <- as.numeric(args[10])
if(is.na(notunegamma)) notunegamma <- 0



print(paste0("np_ratio  is ", np_ratio))
print(paste0("sd beta  is ", sd_dist_beta))
print(paste0("mean beta  is ", mean_dist_beta))
print(paste0("alpha is ", alpha))
print(paste0("sparsity is ", sparsity))
print(paste0("ol percent is ", overlap_pct))
print(paste0("nrep is ", nrep))
print(paste0("p is ", p)) 

# define baseline parameters
n = round(np_ratio * p)

myseed = 2022
set.seed(myseed)

rho_values = c(0.3)
offset = 1
Q = 7
Q_leaves = 4
ending = "hy"


binary_calc_method = "sparse"

# simulation setup - set your paths
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

# LOAD ALL FUNCTIONS 
source(paste0(mydir, "utils/load_aux.R"))

# tell JuliaCall where is the Julia executable (the path can be found by typing `which julia` on terminal)
julia <- julia_setup(JULIA_HOME = my_julia_path)

# load Knockoffs.jl package
julia_library("Knockoffs")
################## GENERATE BETA AND SIGMA #################


########## 1. Generate Group Structure ########

M = 1
group_struct_res <- generate_group_structure(group_size_with_individual = c(1), p, M = 1) 

# groups is a data frame of dimension p x M 
# first column has group number in individual layer 
# second column has group number in second layer 
# ... and so on
groups <- group_struct_res[[1]]

# group_candidates_list is a list whith length == total number of groups (i.e. all individuals + all groups in other layers) 
# each list element contains the number of the individual elements that are part of each group
group_candidates_list <- group_struct_res[[2]]

group_data_frame <- group_struct_res[[3]]


# generate group-size as a list 
group_size_list <- list()

for(m in 1:M) {
  group_size_list[[m]] <- group_data_frame %>% filter(Resolution == m) %>% 
    group_by(SNP) %>%
    count() %>%
    arrange(SNP) %>%
    pull(n)
}

######### 2. Generate Beta ############### 

# Q here denotes the number of leafes!
beta_res <- generate_beta_hierarchical_with_leaf_overlap(p, 
                                                         sparsity, 
                                                         overlap_pct = rep(overlap_pct, 2), 
                                                         Q = Q_leaves, 
                                                         overlap_groups = list(c(1, 2), c(3, 4)), 
                                                         myseed = myseed)

# this is a matrix with 0/1 indicators for whether beta is nonzero
beta_ind_nonzero <- beta_res$beta_ind_nonzero

# non-null for outcome (ANY for global, otherwise must exceed partial conjunction u)
non_null_any <- (rowSums(beta_ind_nonzero) >= 1)

# generate nonzero beta for all global nulls as well 
beta_all_nodes_ind_nonzero <- cbind(beta_ind_nonzero, 
                                    ifelse(rowSums(beta_ind_nonzero[, 1:2]) > 0, 1, 0), 
                                    ifelse(rowSums(beta_ind_nonzero[, 3:4]) > 0, 1, 0), 
                                    ifelse(rowSums(beta_ind_nonzero) > 0, 1, 0))

beta_all_nodes <- beta_all_nodes_ind_nonzero

nonzero_snps_per_q <- list()
nonzero_snps_per_q_character <- c()
for(q in 1:Q) {
  nonzero_snps_per_q[[q]] <- which(beta_all_nodes_ind_nonzero[, q] > 0)
  nonzero_snps_per_q_character <- c(nonzero_snps_per_q_character, paste0("outcome_", q, "_SNP_", nonzero_snps_per_q[[q]]))
}


nonzero_individual_snps <- c()
for(q in 1:Q_leaves) {
  nonzero_individual_snps <- c(nonzero_individual_snps, paste0("outcome_", q, "_SNP_", nonzero_snps_per_q[[q]]))
}


# for each level, calculate fraction for number of outcomes implicated 
fraction_outcome_5 = length(c(nonzero_snps_per_q[[1]], nonzero_snps_per_q[[2]])) / length(unique(c(nonzero_snps_per_q[[1]], nonzero_snps_per_q[[2]])))
fraction_outcome_6 = length(c(nonzero_snps_per_q[[3]], nonzero_snps_per_q[[4]])) / length(unique(c(nonzero_snps_per_q[[3]], nonzero_snps_per_q[[4]])))
fraction_outcome_7 = length(c(nonzero_snps_per_q[[1]], nonzero_snps_per_q[[2]],nonzero_snps_per_q[[3]], nonzero_snps_per_q[[4]])) / length(unique(c(nonzero_snps_per_q[[1]], nonzero_snps_per_q[[2]],nonzero_snps_per_q[[3]], nonzero_snps_per_q[[4]])))

# set gamma
best_gamma = alpha / 2 

################ RUN SIMULATION ##################


all_res <- c()
pct_by_out <- c()
by_resolution <- c()
snr_res <- c()

for(seedA in 1:nrep){
  
  cat(sprintf("Running the %d-th rep.\n", seedA)) 
  
  seedB <- seedA + myseed
  set.seed(seedB)
  
  Y = matrix(NA, nrow = n, ncol = Q)
  
  # generate data for leaves!
  
  
  data_res <- generate_data(Sigma_Dirichlet = FALSE, 
                            distribution_beta = TRUE, 
                            mean_dist_beta = mean_dist_beta,
                            sd_dist_beta = sd_dist_beta,
                            binary_y = TRUE,
                            binary_calc_method = "sparse",
                            beta = beta_ind_nonzero,
                            n = n,
                            p = p,
                            Q = Q_leaves, 
                            rho = rho_values,
                            myseed = seedB, 
                            block_sigma = FALSE, 
                            proportion_of_ones_desired=proportion_of_ones_desired)
  
  
  Y_leaves <- as.matrix(data_res$Y)
  
  Y[, 1:Q_leaves] <- Y_leaves
  Y[, Q_leaves + 1] <- ifelse(rowSums(Y_leaves[, 1:2]) >= 1, 1, 0)
  Y[, Q_leaves + 2] <- ifelse(rowSums(Y_leaves[, 3:4]) >= 1, 1, 0)
  Y[, Q_leaves + 3] <- ifelse(rowSums(Y_leaves) >= 1, 1, 0)
  
  
  X <- as.matrix(data_res$X)
  Sigma <- data_res$Sigma
  
  X.knockoffs <- generate_ko_with_julia(X, Sigma, mu = rep(0, p), method = "maxent", m = 1)$Xko
  
  
  # compute knockoff statistics at each layer for each outcome
  W <- list()
  
  for(q in 1:Q) {
    library(knockoff)
    W[[q]] = knockoff::stat.glmnet_coefdiff(X, X.knockoffs, Y[, q], family = "binomial", cores = 1)
    
  }
  
  
  
  ################## Knockoff rejections: Each Y Separately ######################
  
  taus <- list()
  knockoff_info_list <- list()
  knockoff_fraction <- list()
  
  for(q in 1:Q) {
    print(summary(W[[q]]))
    taus[[q]] = knockoff.threshold(W[[q]], fdr=alpha, offset=offset)
    print(taus[[q]])
    knockoff_info_list[[q]] <- list() 
    knockoff_info_list[[q]]$rejected <- ifelse(W[[q]] > taus[[q]], 1, 0)
    print(sum(knockoff_info_list[[q]]$rejected))
    knockoff_info_list[[q]]$outcome <- rep(q, p)
    knockoff_info_list[[q]]$SNP <- seq(1, p)
    
    knockoff_fraction[[q]] <- (W[[q]] >= taus[[q]]) / (1 + sum(W[[q]] <= -taus[[q]]))
    knockoff_info_list[[q]]$fraction <- knockoff_fraction[[q]]
  }
  
  
  knockoff_rejected <- do.call(rbind.data.frame, knockoff_info_list)
  knockoff_rejected <- knockoff_rejected %>% 
    mutate(subtree = ifelse(outcome %in% c(5, 1, 2), 2, 
                            ifelse(outcome %in% c(3, 4, 6), 3, 1))) %>% 
    mutate(level = ifelse(outcome %in% c(5, 6), (1/2), 
                          ifelse(outcome %in% c(1, 2, 3, 4), 1, (1/4)))) %>%
    mutate(n_outcomes = level) %>%
    mutate(resolution = ifelse(outcome %in% c(5, 6), 2, 
                               ifelse(outcome %in% c(1, 2, 3, 4), 1, 3))) %>%
    filter(rejected == 1)
  
  
  
  
  
  # these knockoff rejections do not have FDR control "by resolution"
  # make sure they have FDR control by resolution! 
  knockoff_info <- do.call(rbind.data.frame, knockoff_info_list)
  knockoff_info <- knockoff_info %>% 
    mutate(subtree = ifelse(outcome %in% c(5, 1, 2), 2, 
                            ifelse(outcome %in% c(3, 4, 6), 3, 1))) %>% 
    mutate(level = ifelse(outcome %in% c(5, 6), (1/2), 
                          ifelse(outcome %in% c(1, 2, 3, 4), 1, (1/4)))) %>%
    mutate(n_outcomes = level) %>%
    mutate(resolution = ifelse(outcome %in% c(5, 6), 2, 
                               ifelse(outcome %in% c(1, 2, 3, 4), 1, 3)))
  
  
  
  if(nrow(knockoff_rejected) > 0) {
    
    knockoff_power_fdp_per_res <- power_fdp_calc_knockoff_per_resolution_hierarchical(knockoff_rejected = knockoff_rejected,
                                                                                      nonzero_snps_per_q_character = nonzero_snps_per_q_character,
                                                                                      nonzero_snps_per_q = nonzero_snps_per_q, 
                                                                                      nonzero_individual_snps = nonzero_individual_snps) 
    
    
    
    knockoff_fdp_by_resolution <- knockoff_power_fdp_per_res$knockoff_fdp_by_resolution
    knockoff_power_by_resolution <- knockoff_power_fdp_per_res$knockoff_power_by_resolution
    knockoff_frac_correct_by_resolution <- knockoff_power_fdp_per_res$knockoff_frac_correct_by_resolution
    knockoff_n_outcomes_implicated = knockoff_power_fdp_per_res$knockoff_n_outcomes_implicated
  } else {
    knockoff_fdp_by_resolution = rep(0, 3)
    knockoff_frac_correct_by_resolution = rep(0, 3)
    knockoff_power_by_resolution = rep(0, 3)
    knockoff_n_outcomes_implicated = rep(0, 3)
  }
  
  
  ################## Knockoff e-BH by resolution ######################
  
  fdp_ebh_knockoff_resolution_all <- c()
  n_outcomes_implicated_ebh_knockoff_resolution_all  <- c()
  frac_correct_rejections_ebh_knockoff_resolution_all <- c()
  power_ebh_knockoff_resolution_all <- c()
  
  for(r in seq(1, 3)) {
    
    
    all_fractions_resolution <- knockoff_info %>% dplyr::filter(resolution == r)
    rejected_fractions_resolution <- knockoff_rejected %>% dplyr::filter(resolution == r)
    
    total_number_groups_in_resolution <- nrow(all_fractions_resolution)
    number_groups_knockoff_rejected <- nrow(rejected_fractions_resolution)
    
    alpha_mod <- alpha * number_groups_knockoff_rejected / total_number_groups_in_resolution
    rejected_fractions_resolution$eval_calc = rejected_fractions_resolution$fraction * (total_number_groups_in_resolution / length(unique(all_fractions_resolution$outcome))) 
    
    ebh_knockoff_resolution_index <- run_eBH(rejected_fractions_resolution$eval_calc, alpha_mod, 
                                             number_groups = number_groups_knockoff_rejected, use_fractions = FALSE)
    
    
    # calculate power / fdr
    # calculate power / fdr after running ebh on knockoff outer nodes
    if(length(ebh_knockoff_resolution_index) > 0 & !is.na(ebh_knockoff_resolution_index[1])) {
      
      ebh_rejected_fractions_resolution = rejected_fractions_resolution[ebh_knockoff_resolution_index, ]
      
      ebh_knockoff_resolution_power_fdp <- power_fdp_calc_hierarchical(nonzero_snps_per_q = nonzero_snps_per_q, 
                                                                       all_outcomes = seq(1, 7) ,
                                                                       rejections = ebh_rejected_fractions_resolution, 
                                                                       nonzero_snps_per_q_character = nonzero_snps_per_q_character, 
                                                                       nonzero_individual_snps = nonzero_individual_snps)
      
      fdp_ebh_knockoff_resolution = ebh_knockoff_resolution_power_fdp$fdp 
      power_ebh_knockoff_resolution = ebh_knockoff_resolution_power_fdp$power
      frac_correct_rejections_ebh_knockoff_resolution = ebh_knockoff_resolution_power_fdp$frac_correct_rejections
      n_outcomes_implicated_ebh_knockoff_resolution = ebh_knockoff_resolution_power_fdp$n_outcomes_implicated
      
    } else {
      fdp_ebh_knockoff_resolution <- 0
      power_ebh_knockoff_resolution <- 0
      frac_correct_rejections_ebh_knockoff_resolution = 0
      n_outcomes_implicated_ebh_knockoff_resolution = 0
    }
    
    
    fdp_ebh_knockoff_resolution_all <- c(fdp_ebh_knockoff_resolution_all, fdp_ebh_knockoff_resolution)
    power_ebh_knockoff_resolution_all <- c(power_ebh_knockoff_resolution_all, power_ebh_knockoff_resolution)
    frac_correct_rejections_ebh_knockoff_resolution_all <- c(frac_correct_rejections_ebh_knockoff_resolution_all, frac_correct_rejections_ebh_knockoff_resolution)
    n_outcomes_implicated_ebh_knockoff_resolution_all <- c(n_outcomes_implicated_ebh_knockoff_resolution_all, n_outcomes_implicated_ebh_knockoff_resolution)
    
    
  }
  
  
  ################## Knockoff outer nodes ######################
  
  # filter to knockoff "outer nodes" for each outcome
  # for each "subtree" or root-node, filter to max outcome 
  # this keeps the mnost specific rejected node among each of the "sub-trees" (1, 2, 5), (3, 4, 6), (7 = all)
  knockoff_rejected_outer <- knockoff_rejected %>% 
    group_by(SNP, subtree) %>% 
    dplyr::filter(level == max(level))
  
  # for each snp, check whether another outcome was rejected that is not Y7
  knockoff_max_subtree_by_snp <- knockoff_rejected_outer %>% 
    group_by(SNP) %>% 
    mutate(max_subtree = max(subtree)) %>% 
    dplyr::select(SNP, max_subtree) %>% 
    unique()
  
  knockoff_rejected_outer <- left_join(knockoff_rejected_outer, knockoff_max_subtree_by_snp, by = "SNP")
  
  # if Y1 was rejected but also a node further down the tree, only keep the node further down
  # however, we can keep rejected nodes in different subtrees
  knockoff_rejected_outer <- knockoff_rejected_outer %>% 
    dplyr::filter(!(subtree == 1 & max_subtree > 1))
  
  if(nrow(knockoff_rejected_outer) > 0) {
    
    knockoff_power_fdp <- power_fdp_calc_hierarchical(nonzero_snps_per_q = nonzero_snps_per_q, 
                                                      all_outcomes = seq(1, 7) ,
                                                      rejections = knockoff_rejected_outer, 
                                                      nonzero_snps_per_q_character = nonzero_snps_per_q_character, 
                                                      nonzero_individual_snps = nonzero_individual_snps)
    
    fdp_knockoff = knockoff_power_fdp$fdp 
    power_knockoff = knockoff_power_fdp$power
    frac_correct_rejections_knockoff = knockoff_power_fdp$frac_correct_rejections
    n_outcomes_implicated_knockoff = knockoff_power_fdp$n_outcomes_implicated
    
    
    print("number rejections knockoff")
    print(nrow(knockoff_rejected_outer))
    
    
  } else {
    fdp_knockoff <- 0
    power_knockoff <- 0
    frac_correct_rejections_knockoff = 0
    n_outcomes_implicated_knockoff = 0
  }
  
  print(fdp_knockoff) 
  
  ################## Knockoff outer nodes + e-BH ######################
  
  total_number_groups <- nrow(knockoff_info)
  number_groups_outer <- nrow(knockoff_rejected_outer)
  alpha_mod <- alpha * number_groups_outer / total_number_groups
  knockoff_rejected_outer$eval_calc = knockoff_rejected_outer$fraction * (total_number_groups / Q) # every fraction gets multiplied by p = 1000
  
  
  ebh_knockoff_outer_index <- run_eBH(knockoff_rejected_outer$eval_calc, alpha_mod, 
                                      number_groups = number_groups_outer, use_fractions = FALSE)
  
  
  # calculate power / fdr after running ebh on knockoff outer nodes
  if(length(ebh_knockoff_outer_index) > 0 & !is.na(ebh_knockoff_outer_index[1])) {
    
    ebh_knockoff_outer_rejections = knockoff_rejected_outer[ebh_knockoff_outer_index, ]
    
    ebh_knockoff_outer_power_fdp <- power_fdp_calc_hierarchical(nonzero_snps_per_q = nonzero_snps_per_q, 
                                                                all_outcomes = seq(1, 7) ,
                                                                rejections = ebh_knockoff_outer_rejections, 
                                                                nonzero_snps_per_q_character = nonzero_snps_per_q_character, 
                                                                nonzero_individual_snps = nonzero_individual_snps)
    
    fdp_ebh_knockoff = ebh_knockoff_outer_power_fdp$fdp 
    power_ebh_knockoff = ebh_knockoff_outer_power_fdp$power
    frac_correct_rejections_ebh_knockoff = ebh_knockoff_outer_power_fdp$frac_correct_rejections
    n_outcomes_implicated_ebh_knockoff = ebh_knockoff_outer_power_fdp$n_outcomes_implicated
    
    
    print("number rejections knockoff")
    print(nrow(ebh_knockoff_outer_rejections))
    
    
  } else {
    fdp_ebh_knockoff <- 0
    power_ebh_knockoff <- 0
    frac_correct_rejections_ebh_knockoff = 0
    n_outcomes_implicated_ebh_knockoff = 0
  }
  
  
  
  #################   Get e-values ############################
  
  taus_evals <- list() 
  fracs <- list()
  evals <- list() 
  
  # get e-values for each outcome separately
  for(q in 1:Q) {
    taus_evals[[q]] = knockoff.threshold(W[[q]], fdr=best_gamma, offset=offset)
    print(taus_evals[[q]])
    
    evals[[q]] <- p * ((W[[q]] >=  taus_evals[[q]] ) / (1 + sum(W[[q]] <= -taus_evals[[q]] )))
    fracs[[q]] <- list() 
    fracs[[q]]$fracs <- ((W[[q]] >=  taus_evals[[q]] ) / (1 + sum(W[[q]] <= -taus_evals[[q]] )))
    fracs[[q]]$outcome <- rep(q, p)
    fracs[[q]]$SNP <- seq(1, p)
  }
  
  
  # now put through focused e-BH procedure
  fracs_df <- do.call(rbind.data.frame, fracs)
  fracs_df <- fracs_df %>% 
    mutate(subtree = ifelse(outcome %in% c(5, 1, 2), 2, 
                            ifelse(outcome %in% c(3, 4, 6), 3, 1))) %>% 
    mutate(level = ifelse(outcome %in% c(5, 6), (1/2), 
                          ifelse(outcome %in% c(1, 2, 3, 4), 1, (1/4)))) %>% 
    mutate(n_outcomes = level)
  
  ################## run e-FBH #######################
  
  cat("running focused eBH")
  focusedeBH_rejections <- focusedeBH_hierarchical_y(fracs = fracs_df$fracs,
                                                     fracs_df = fracs_df,
                                                     numerator_threshold = Q,
                                                     alpha = 0.2) 
  cat("done\n")
  
  
  if(!is.null(focusedeBH_rejections)) {
    # calculate power and fdr 
    focusedeBH_power_fdp <- power_fdp_calc_hierarchical(nonzero_snps_per_q = nonzero_snps_per_q,
                                                        all_outcomes = seq(1, 7), 
                                                        rejections = focusedeBH_rejections, 
                                                        nonzero_snps_per_q_character = nonzero_snps_per_q_character, 
                                                        nonzero_individual_snps = nonzero_individual_snps)
    
    fdp_fbh = focusedeBH_power_fdp$fdp 
    power_fbh = focusedeBH_power_fdp$power
    frac_correct_rejections_fbh = focusedeBH_power_fdp$frac_correct_rejections
    n_outcomes_implicated_fbh = focusedeBH_power_fdp$n_outcomes_implicated
    
    print("fdp, power, frac")
    print(fdp_fbh)
    print(power_fbh)
    print(frac_correct_rejections_fbh)
    
    pct_rejected_by_outcome <- c()
    n_rejected = nrow(focusedeBH_rejections)
    for(q in 1:Q) {
      n_rejected_q <- nrow(focusedeBH_rejections %>% dplyr::filter(outcome == q))
      pct_rejected_by_outcome[q] <- n_rejected_q / n_rejected
    }
    
  } else {
    print("no rejections")
    power_fbh = 0
    fdp_fbh = 0
    frac_correct_rejections_fbh = 0
    n_outcomes_implicated_fbh = 0
    pct_rejected_by_outcome <- rep(NA, Q) # this is given there were some rejections 
    # do not want to count as 0
  }
  
  
  # combine output
  summary_df <- data.frame(method = c("knockoffouter", "ebhknockoffouter","eFBH"), 
                           power = c(power_knockoff, power_ebh_knockoff, power_fbh), 
                           fdr = c(fdp_knockoff, fdp_ebh_knockoff, fdp_fbh), 
                           frac_correct_rejections = c(frac_correct_rejections_knockoff, frac_correct_rejections_ebh_knockoff, frac_correct_rejections_fbh),
                           n_outcomes_implicated = c(n_outcomes_implicated_knockoff,n_outcomes_implicated_ebh_knockoff, n_outcomes_implicated_fbh),
                           np_ratio = np_ratio, 
                           seed = seedB, 
                           alpha = alpha, 
                           overlap_pct = overlap_pct,
                           gamma = best_gamma)
  
  
  summary_pct_rejected = data.frame(outcome = paste0("Y", seq(1, Q)), 
                                    pct_rejected = pct_rejected_by_outcome, 
                                    np_ratio = np_ratio, 
                                    seed = seedB, 
                                    alpha = alpha, 
                                    overlap_pct = overlap_pct,
                                    gamma = best_gamma)
  
  
  summary_by_resulution = data.frame(resolution = rep(seq(1, 3), 2), 
                                     method = c(rep("knockoff", 3), rep("ebh knockoff", 3)),
                                     fdp = c(knockoff_fdp_by_resolution, fdp_ebh_knockoff_resolution_all), 
                                     power = c(knockoff_power_by_resolution, power_ebh_knockoff_resolution_all), 
                                     frac_correct = c(knockoff_frac_correct_by_resolution, frac_correct_rejections_ebh_knockoff_resolution_all), 
                                     n_outcomes_implicated = c(knockoff_n_outcomes_implicated, n_outcomes_implicated_ebh_knockoff_resolution_all), 
                                     np_ratio = np_ratio, 
                                     seed = seedB, 
                                     alpha = alpha, 
                                     overlap_pct = overlap_pct, 
                                     gamma = best_gamma)
  
  all_res <- rbind(all_res, summary_df)
  pct_by_out <- rbind(pct_by_out, summary_pct_rejected)
  by_resolution <- rbind(by_resolution, summary_by_resulution)
  
  
}  


all_res %>% group_by(method) %>% summarise(mean_frac = mean(frac_correct_rejections))
by_resolution %>% group_by(method, resolution) %>% summarise(mean_frac = mean(frac_correct))


out_dir <- sprintf("%s/hp_union_by_res_%.1f_ol_%.2f_sp_%.4f_nrep_%d_r%s_%s_al%.4f_dpo%.2f_sdb%.2f_mdb%.2f_p%d_ntg%d.csv", save_dir, 
                   np_ratio, overlap_pct, sparsity, nrep, paste(rho_values, collapse = ""), ending, 
                   alpha, proportion_of_ones_desired, sd_dist_beta, mean_dist_beta, p, notunegamma)
write_csv(data.frame(by_resolution), out_dir)

out_dir <- sprintf("%s/hp_union_%.1f_ol_%.2f_sp_%.4f_nrep_%d_r%s_%s_al%.4f_dpo%.2f_sdb%.2f_mdb%.2f_p%d_ntg%d.csv", save_dir, 
                   np_ratio, overlap_pct, sparsity, nrep, paste(rho_values, collapse = ""), ending, 
                   alpha, proportion_of_ones_desired, sd_dist_beta, mean_dist_beta, p, notunegamma)
write_csv(data.frame(all_res), out_dir)

out_dir <- sprintf("%s/hp_union_by_out_%.1f_ol_%.2f_sp_%.4f_nrep_%d_r%s_%s_al%.4f_dpo%.2f_sdb%.2f_mdb%.2f_p%d_ntg%d.csv", save_dir, 
                   np_ratio, overlap_pct, sparsity, nrep, paste(rho_values, collapse = ""), ending,
                   alpha, proportion_of_ones_desired, sd_dist_beta, mean_dist_beta, p, notunegamma)
write_csv(data.frame(pct_by_out), out_dir)

