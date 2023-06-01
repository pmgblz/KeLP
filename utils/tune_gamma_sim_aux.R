

                # Tune gamma in simulations#


# function to tune gamma in the different scenarios: 
  # Kelp / Focused e-BH 
  # Hierarchical Outcomes 
  # E-MLKF



# function to tune gamma on independently generated training set with the same specifications 
# for each gamma in gamma sequence, repeat data generation process gamma_reps time and calculate the number of rejections with eblipr 
# average the number of eblipr rejections across all gamma_reps repetitions
# pick the gamma that leads to the largest number of rejections with eblipr
# INPUT: 
  # gamma_seq: vector of gammas to try
  # gamma reps: How often to generate the training data set. Results will be averaged across these datasets. 
  # n: number of observations. 
  # p: number of predictors.
  # amp: signal amplitude.
  # Sigma: covariance matrix. 
  # beta: matrix of betas, dimension n x Q
  # Q: number of outputs
  # groups: data frame of dimension p x M containing group information
  # group_candidates_list: list of of length equal to total number of groups, each list element contains information about the group 
  # M: number of layers
  # statistic_type: type of importance statistic to be calculated
  # group_size_with_individual: vector containing group size, including individual layer
  # myseed: seed
  # use_amp: generate data using amplitude or SNR
  # knockoff_method: method to construct knockoffs 

# OUTPUT: 
  # best_gamma: gamma resulting in the largest number of rejections






tune_gamma <- function(gammas, 
                       groups_combined, 
                       group_size_with_individual,
                       Sigma_Dirichlet = FALSE, 
                       distribution_beta = FALSE, 
                       binary_y = FALSE,
                       intercept_binary = FALSE,
                       beta,
                       n,
                       p,
                       Q, 
                       k, 
                       alpha0_dirichlet, 
                       alphasum_dirichlet, 
                       max_corr, 
                       rho,
                       amp,
                       mean_dist_beta,
                       sd_dist_beta,
                       myseed, 
                       block_sigma, 
                       block_size_vector) {
  
  
  set.seed(myseed)
  
  # GENERATE DATA ONCE 
  data_res <- generate_data(Sigma_Dirichlet = Sigma_Dirichlet, 
                            block_sigma = block_sigma,
                            distribution_beta = distribution_beta, 
                            amp = amp,
                            binary_y = binary_y,
                            beta = beta,
                            n = n,
                            p = p,
                            rho = rho,
                            Q = Q, 
                            k = k, 
                            alpha0_dirichlet = alpha0_dirichlet, 
                            alphasum_dirichlet = alphasum_dirichlet, 
                            max_corr = max_corr, 
                            mean_dist_beta = mean_dist_beta,
                            sd_dist_beta = sd_dist_beta,
                            myseed = myseed, 
                            groups_combined, 
                            group_size_with_individual = group_size_with_individual, 
                            block_size_vector = block_size_vector)
  
  print(summary(data_res$Y))
 
  
  Y <- data_res$Y
  X <- data_res$X
  
  etas <- data_res$etas
  Sigma <- data_res$Sigma
  
  
  # define objects to hold results
  fracs_list <- list()
  nonzero_group_per_layer <- list()
  nonzero_res_group_per_layer <-list()

  cand_groups <- list()
  
  counter_layer = 1
  
  
  W_list <- list()
  
  ############### ... BEGIN CUTOFF LOOP #############
  for(c in 1:M) {
    
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
    
    
    
    #for(gamma in c(alpha / 2)) { alpha / 20, alpha / 15, alpha / 10, alpha / 6, alpha / 4, alpha / 2, 
    
    
    fracs_per_q <- matrix(NA, nrow = length(unique(groups)), ncol = Q)
    selected_knockoff_groups_per_q <- list()
    
    
    ############### ...... BEGIN OUTCOME LOOP #############
    
    for(q in 1:Q) {
      
      y = as.matrix(Y[, q])
      
      if((c == 1) | (number_unique_groups == p)) {
        result <- run_ko_lasso_with_julia(y = y, X = X, Sigma = Sigma, mu = mu, method = method, m = 1, alpha = alpha)
        
      } else {
        # call Julia solver
        result <- run_group_ko_lasso_with_julia(y = y, X = X, Sigma = Sigma, mu = mu, method = method, m = 1, groups = groups, alpha = alpha)
      }
      
      
      
      # get rejected groups by the knockoff 
      selected_knockoff_groups_per_q[[q]] <- result$selected[1]
      
      # get fractions for e-values
      W_list[[c]][[q]] <- result$W
      
      
    }
    
    # increase counter for next layer
    counter_layer = counter_layer + 1
  }
  
  
  
  # THEN LOOP ACROSS GAMMAS, SAME SEED
  
  all_res_gammas <- c()
  all_res_snps <- c()
  opt_value_combined <- c()
  number_rejections <- c()
  for(r in 1:nrow(gammas)) {
    
    
    
    
    # begin loop to construct e-values
    
    for(c in 1:M) {
      
      fracs_per_q <- matrix(NA, nrow = length(unique(groups_combined[, c(c)])), ncol = Q)
      
      for(q in 1:Q) {
        
        tau_gamma <- knockoff.threshold(W_list[[c]][[q]], fdr=gammas[r, c], offset=1)
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
    
    # get eblipr rejections    
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
      cat("Running eblipr ...")
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
      eblipr_detections <- res[[1]]
      
      # get selected groups
      selected_eblipr_groups <- sapply(eblipr_detections, `[[`, "Res_Group")
      
      selected_eblipr_fracs <- sapply(eblipr_detections, `[[`, "fracs")
      included_eblipr_SNPs <- unlist(sapply(eblipr_detections, `[[`, "SNP"))
      
    } else {
      selected_eblipr_groups <- NULL
      selected_eblipr_fracs <- NULL
      included_eblipr_SNPs <- NULL
    }
    
    # get power and fdr for kelp  # check if self consistent
    
    if(length(selected_eblipr_groups) > 0) {
      # get frac correct
      eblipr_frac_correct = length(intersect(included_eblipr_SNPs, nonzero_any_indiv)) / length(nonzero_any_indiv)
    } else {
      eblipr_frac_correct = 0
    }
    
    # 
    if(length(selected_eblipr_groups) > 0) {
      opt_value <- eblipr_frac_correct / (length(included_eblipr_SNPs))
    } else {
      opt_value = eblipr_frac_correct 
    }
    
    opt_value_combined <- c(opt_value_combined, opt_value )
    all_res_snps <- c(all_res_snps, length(included_eblipr_SNPs))
    all_res_gammas <- c(all_res_gammas, eblipr_frac_correct)
    number_rejections <- c(number_rejections, length(selected_eblipr_groups))
    
  
    
    
  }
  
  # 
  # # 0.1 = alpha / 2
  # if(0.1 %in% gammas[all_res_gammas == max(all_res_gammas), ]) {
  #   max_gamma = rep(0.1, )
  # } else {
    #max_gamma <- gammas[tail(which(all_res_gammas==max(all_res_gammas), arr.ind=TRUE) , 1), ]
  #max_gamma <- gammas[which.max(opt_value_combined), ]
  # }

 

  return(list(all_res_gammas = all_res_gammas, 
              all_res_snps = all_res_snps, 
              opt_value_combined = opt_value_combined, 
              number_rejections = number_rejections))
  
}



# tune hierarchical gamma 

tune_gamma_hierarchical <- function(gammas_tune, 
                                    beta_all_nodes, 
                                    n, 
                                    p, 
                                    Q, 
                                    rho_values, 
                                    prevalance_log,
                                    sd_dist_beta,
                                    mean_dist_beta,
                                    myseed, 
                                    nonzero_individual_snps = nonzero_individual_snps, 
                                    union = FALSE, 
                                    Q_leaves, 
                                    binary_calc_method = binary_calc_method, 
                                    proportion_of_ones_desired, 
                                    fraction_outcome_5, 
                                    fraction_outcome_6, 
                                    fraction_outcome_7) {
  
  
  
  # generate data
  if(!union) {
    data_res_tune <- generate_data(Sigma_Dirichlet = FALSE, 
                                   distribution_beta = TRUE, 
                                   sd_dist_beta = sd_dist_beta, 
                                   mean_dist_beta = mean_dist_beta,
                                   binary_y = TRUE,
                                   binary_calc_method = "default",
                                   block_sigma = FALSE,
                                   beta = beta_all_nodes,
                                   n = n,
                                   p = p,
                                   Q = Q, 
                                   rho = rho_values,
                                   amp = amp,
                                   nonzero_snps_per_q, 
                                   nonzero_snps_per_q_character, 
                                   myseed = myseed)
    
    
    
    Y_tune <- as.matrix(data_res_tune$Y)
    
  } 
  
  if(union) {
    data_res_tune <- generate_data(Sigma_Dirichlet = FALSE, 
                                   distribution_beta = TRUE, 
                                   sd_dist_beta = sd_dist_beta, 
                                   mean_dist_beta = mean_dist_beta,
                                   binary_y = TRUE,
                                   binary_calc_method = "sparse",
                                   block_sigma = FALSE,
                                   beta = beta_all_nodes,
                                   n = n,
                                   p = p,
                                   Q = Q_leaves, 
                                   rho = rho_values,
                                   binary_y_intercept = -log(prevalence_log),
                                   amp = amp,
                                   nonzero_snps_per_q, 
                                   nonzero_snps_per_q_character, 
                                   myseed = myseed, 
                                   proportion_of_ones_desired = proportion_of_ones_desired)
    
    
    Y_tune <- matrix(NA, nrow = n, ncol = Q)
    
    Y_tune_leaves <- as.matrix(data_res_tune$Y)
    Y_tune[, 1:Q_leaves] <- Y_tune_leaves
    Y_tune[, Q_leaves + 1] <- ifelse(rowSums(Y_tune_leaves[, 1:2]) >= 1, 1, 0)
    Y_tune[, Q_leaves + 2] <- ifelse(rowSums(Y_tune_leaves[, 3:4]) >= 1, 1, 0)
    Y_tune[, Q_leaves + 3] <- ifelse(rowSums(Y_tune_leaves) >= 1, 1, 0)
  }
  
  X_tune <- as.matrix(data_res_tune$X)
  Sigma <- data_res_tune$Sigma
  
  #diags <- knockoff::create.solve_asdp(Sigma)
  #X.knockoffs.tune <- create.gaussian(X_tune, rep(0,p), Sigma, diag_s = diags)
  X.knockoffs.tune <- generate_ko_with_julia(X_tune, Sigma, mu = rep(0, p), method = "maxent", m = 1)$Xko
  
  
  # compute knockoff statistics at each layer for each outcome
  W.tune <- list()
  
  for(q in 1:Q) {
    library(knockoff)
    W.tune[[q]] = knockoff::stat.glmnet_coefdiff(X_tune, X.knockoffs.tune, Y_tune[, q], family = "binomial", cores = 1)
    
  }
  
  for(i in 1:length(gammas_tune)){
    
    # get evals
    taus_evals <- list()
    fracs <- list()
    evals <- list()
    
    # get e-values for each outcome separately
    for(q in 1:Q) {
      taus_evals[[q]] = knockoff.threshold(W.tune[[q]], fdr=gammas_tune[i], offset=1)
      print(taus_evals[[q]])
      
      evals[[q]] <- p * ((W.tune[[q]] >=  taus_evals[[q]] ) / (1 + sum(W.tune[[q]] <= -taus_evals[[q]] )))
      fracs[[q]] <- list()
      fracs[[q]]$fracs <- ((W.tune[[q]] >=  taus_evals[[q]] ) / (1 + sum(W.tune[[q]] <= -taus_evals[[q]] )))
      fracs[[q]]$outcome <- rep(q, p)
      fracs[[q]]$SNP <- seq(1, p)
    }
    
    
    # now put through focused e-BH procedure
    fracs_df <- do.call(rbind.data.frame, fracs)
    fracs_df <- fracs_df %>%
      mutate(subtree = ifelse(outcome %in% c(5, 1, 2), 2,
                              ifelse(outcome %in% c(3, 4, 6), 3, 1))) %>%
      mutate(n_outcomes = ifelse(outcome %in% c(1, 2, 3, 4), 1, 
                                        ifelse(outcome == 5, 2*fraction_outcome_5, 
                                               ifelse(outcome == 6, 2*fraction_outcome_6, 
                                                      4*fraction_outcome_7)))) %>%
      mutate(level = ifelse(outcome %in% c(5, 6), (1/2),
                            ifelse(outcome %in% c(1, 2, 3, 4), 1, (1/4))))
    
    
    
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
      
      power_fbh = focusedeBH_power_fdp$power
      
      
    } else {
      power_fbh = 0
    }
    
    all_gamma_tune_res <- c(all_gamma_tune_res, power_fbh)
    
    
  }
  
  # choose 0.1 if among all that are equal
  if((alpha / 2) %in% gammas_tune[all_gamma_tune_res == max(all_gamma_tune_res)]) {
    best_gamma = alpha / 2
  } else {
    best_gamma <- gammas_tune[which.max(all_gamma_tune_res)]
  }
  
  
  return(best_gamma)
}


# tune gamma for MLKF 
# parts of code-structure taken from Katsevich, Sabatti (2019) Multilayer knockoff filter: Controlled variable selection at multiple resolutions; downloaded from https://katsevich-lab.github.io/publications/
# if gamma_m seq params would create negative alpha_emlkf, starts at 0 (same as min(1, gamma_m max))

# INTPUT: 


  # X
  # Y 
  # W: importance statistics; list of length M 
  # groups: data frame of dimension p x M containing group information 
  # alpha_emlkf: target FDR levels 
  # gamma_m: starting gamma used for each layer 
  # FDP_hat_type: regulates offset 
  # gamma_seq_params: parameters automatically generating gamma sequence based on gamma_m 
  # same_gamma_each_m: indicator whether the same gamma should be used in every layer
  # tail_max: Indicator for whether to return the largest or smallest index when ties 
emlkf_tune <-  function(X, Y, W, 
                        groups, 
                        alpha_emlkf, 
                        gamma_m, 
                        FDP_hat_type, 
                        gamma_seq_params = c(0.1, 0.2, 0.001), 
                        same_gamma_each_m = FALSE, 
                        tail_max = FALSE, 
                        myseed){
  
  
  Sys.sleep(1)
  set.seed(myseed)
  Sys.sleep(1)
  
  # problem dimensions
  N = nrow(X) # number of samples
  n = nrow(groups) # number of variables
  M = ncol(groups) # number of layers
  G = apply(groups,2,max) # G[m] = number of groups for at layer m
  
  # check input for correctness
  stopifnot(length(Y) == N)
  stopifnot(ncol(X) == n)
  stopifnot(length(alpha_emlkf) == M)
  stopifnot(FDP_hat_type %in% c("kn", "kn+"))
  
  if(FDP_hat_type == "kn") {
    my_offset = 0
  }
  
  if(FDP_hat_type == "kn+") {
    my_offset = 1
  }
  
  
  # options find best gamma jointly for all m by either 
  #               (i) forcing gamma to be the same in each layer or 
  #               (ii) considering all "combinations" of the gamma-sequence 
  
  
  #  one gamma for all layers
  if(same_gamma_each_m) {
    
    
    # assumes all gammas are the same
    combs_seq =  unique(c(seq(max(0, gamma_m[[1]] - gamma_seq_params[1]), gamma_m[[1]], by = gamma_seq_params[3]),
                          seq(gamma_m[[1]], min(1, gamma_m + gamma_seq_params[2]), by = gamma_seq_params[3])))
    
    combs <- data.frame()
    
    for(s in 1:length(combs_seq)) {
      combs <- rbind(combs, rep(combs_seq[s], M))
    }
    colnames(combs) <- paste0("Layer",seq(1:M))
    
  } else {  # different gamma for each m 
    
    # get all alpha_emlkf combinations based on broad sequence
    combs <- list()
    
    for(m in 1:M) {
      combs[[m]] <- unique(c(seq(max(0, gamma_m[m] - gamma_seq_params[1]), gamma_m[m], by = gamma_seq_params[3]),
                             seq(gamma_m[m], min(1, gamma_m[m] + gamma_seq_params[2]), by = gamma_seq_params[3])))
    }
    
    names(combs) <- paste0("Layer",seq(1:M))
    combs <- combs %>% cross_df()
    
    
  }
  
  print(combs)
  cat("Find best gamma combination ...")
  res_gammas_search <- apply(combs, 1, find_num_rejections,
                             alpha_emlkf = alpha_emlkf, 
                             groups = groups, 
                             my_offset = my_offset, 
                             G = G, W = W, 
                             myseed = myseed) # returns a list
  cat("done.\n")
  
  # convert the output of res_gammas_search to df
  all_num_rejections <- c()
  taus_gamma <- matrix(NA, nrow = nrow(combs), ncol = M)
  for(i in 1:nrow(combs)){
    all_num_rejections[i] <- res_gammas_search[[i]]$num_rejections
    taus_gamma[i, ] <- res_gammas_search[[i]]$taus_gamma
  }
  
  
  # get the combination which leads to largest number of rejections
  if(tail_max) {
    final_gamma_combination <- as.matrix(combs[tail(which(all_num_rejections==max(all_num_rejections)),1), ])
  } else {
    
    final_gamma_combination <- as.matrix(combs[which.max(all_num_rejections), ])
    
    all_comb_max_rejections <- as.matrix(combs[which(all_num_rejections==max(all_num_rejections)), ])
    
    # if for all ties it contains alpha / 2, choose this 
    print(all_comb_max_rejections)
    if((alpha / 2) %in% all_comb_max_rejections[, 1] | (alpha / 2) %in% all_comb_max_rejections[, 2]) {
      
      # if alpha / 2 is in both l1 and l2, then choose the column that has most
      if((alpha / 2) %in% all_comb_max_rejections[, 2] & ((alpha / 2) %in% all_comb_max_rejections[, 1])) {
        final_gamma_combination <- all_comb_max_rejections[ all_comb_max_rejections[, 1] == (alpha / 2) & all_comb_max_rejections[, 2] == (alpha / 2), ]
      } 
      
      if((alpha / 2) %in% all_comb_max_rejections[, 1] & !( (alpha / 2) %in% all_comb_max_rejections[, 2])) {
        final_gamma_combination <- all_comb_max_rejections[ all_comb_max_rejections[, 1] == (alpha / 2), ]
        # if there are multiple chooose the max
        if(!is.null(dim(final_gamma_combination)[1])) {
          final_gamma_combination <- all_comb_max_rejections[ all_comb_max_rejections[, 1] == (alpha / 2), ][dim(final_gamma_combination)[1], ]
        }
      } 
      
      if((alpha / 2) %in% all_comb_max_rejections[, 2] & !( (alpha / 2) %in% all_comb_max_rejections[, 1])) {
        final_gamma_combination <- all_comb_max_rejections[ all_comb_max_rejections[, 2] == (alpha / 2), ]
        # if there are multiple chooose the max
        if(!is.null(dim(final_gamma_combination)[1])) {
          final_gamma_combination <- all_comb_max_rejections[ all_comb_max_rejections[, 2] == (alpha / 2), ][dim(final_gamma_combination)[1], ]
        }
      }
     
    } else {
      final_gamma_combination <- as.matrix(combs[which.max(all_num_rejections), ])
    }
    
  }
  
 
  print(all_num_rejections)
  print(max(all_num_rejections))

  print(final_gamma_combination)
  
  return(list(final_gamma_combination = final_gamma_combination, 
              all_num_rejections = all_num_rejections, 
              combs = combs))
  
}


emlkf_return_gamma <- function(X, Y, n, p, amp, Sigma, beta, Q, groups, knockoff_type,  alpha_emlkf, 
                               starting_gamma, FDP_hat_type, gamma_seq_params = c(0.1, 0.1, 0.005), 
                               same_gamma_each_m = TRUE, tail_max = TRUE, 
                               myseed, 
                               generate_data = FALSE, rho) {
  
    print("TUNING GAMMA")
    
    # note that this is a different seed than the one used in the simulations to ensure that "noise" is different
  Sys.sleep(1)
  set.seed(myseed)
  Sys.sleep(1)
    
    
    if(generate_data) {
      # generate data
      data_res <- generate_data(Sigma_Dirichlet = FALSE, 
                                block_sigma = FALSE,
                                distribution_beta = FALSE, 
                                binary_y = FALSE,
                                beta = beta,
                                Q = Q, 
                                n = n,
                                p = p, 
                                rho = rho_values, 
                                amp = amp, 
                                myseed = myseed)
      Y <- as.matrix(data_res$Y)
      X <- as.matrix(data_res$X)
      Sigma <- as.matrix(data_res$Sigma)
      
      print(summary(Y))
      print(summary(rowMeans(X)))
    }
    
    
    # standardize columns of X, as in MLKF
    X = scale(X, center = FALSE)/sqrt(n-1)

    # construct knockoffs at each layer without Julia
    Sys.sleep(1)
    set.seed(myseed)
    Sys.sleep(1)
    X.knockoffs = get_knockoff_variables(X, groups, knockoff_type)
    
    print("knockoff summary")
    print(summary(rowMeans(X.knockoffs[[1]])))
    print(summary(rowMeans(X.knockoffs[[2]])))
    
    # compute knockoff statistics at each layer for each outcome
    W <- list()
    
    for(q in 1:Q) {
      W[[q]] = get_knockoff_stats(X, X.knockoffs, Y[, q], groups, statistic_type)
    }
    
    
    print("W summary")
    print(summary((W[[1]][[1]])))
    print(summary((W[[1]][[2]])))
    
    
    # find best gamma (choose the one that leads to the largest number of rejections)
    # same gamma in each layer: in simulations it looks like when allowing gamma to be different in each layer, there are 
    # repetitive combinations that allow for the largest number of rejections, intermixed with zeros 
    # therefore, it seems to be possible to just have gamma the same in each layer
    output_e_mlkf <- emlkf_tune(X = X, Y = Y, W = W[[Q]], groups = groups, 
                                alpha_emlkf = alpha_emlkf,
                                gamma_m = rep(starting_gamma, M), 
                                FDP_hat_type = FDP_hat_type, 
                                gamma_seq_params = gamma_seq_params, #c(0.1, 0.1, 0.001) 
                                same_gamma_each_m = same_gamma_each_m, 
                                tail_max = tail_max, 
                                myseed = myseed)
    
    final_gamma_combination <- output_e_mlkf$final_gamma_combination
    print("FINAL GAMMA COMBINATION IS:")
    print(final_gamma_combination)
    
    return(final_gamma_combination)

  
}


