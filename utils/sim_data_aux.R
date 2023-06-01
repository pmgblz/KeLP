

# Data Generation Functions #


# This file contains functions to simulate data. 
# Specifically this file contains functions for the following three parts: 


# 1. Generate X and Y matrices based on Beta and Sigma
# 2. Generate Beta
# 3. Generate Sigma
# 4. Generate desired group structure
# 5. Generate hierarchically structured outputs


########################## GENERATE GROUP STRUCTURE #######################


# function to generate group structure

# INPUT: 
# group_size_with_individual: vector containing the size of the groups, including individual layer with size = 1 
# p: number of predictors 
# M: number of layers
# OUTPUT: group membership information in the following formats:
# groups: data frame of dimensions p x M specifying group number for every individual 
# group_candidates_list: list containing group information. each list element is a group-list, 
# with elements SNP (containing the snps or individual hypotheses involved), 
# the group number, the Resolution, and a Resolution-group identifier
# groups_data_frame: data frame containing the information in group_candidates_list in data-frame format
generate_group_structure <- function(group_size_with_individual, p, M) {
  
  
  # create two layers of groups, one singleton layer, and one layer containing 25 groups of size 10 each
  groups = matrix(0, p, M)
  groups[,1] = 1:p
  if(M > 1) {
    for(m in 2:M) {
      groups[,m] = rep(1:(p/group_size_with_individual[m]), each = group_size_with_individual[m])
    }
  }
  
  # based on group matrix, construct candidate groups
  group_candidates_list <- list()
  for(m in 1:M) {
    temp <- split(groups[, 1],groups[, m])
    group_number <- paste0("group_", names(temp))
    temp <- unname(temp)
    Resolution = lapply(temp,
                        function(x) which(group_size_with_individual == length(x)))
    
    group_candidates_list <- c(group_candidates_list, mapply(list, SNP = temp, group = group_number, Resolution = Resolution, SIMPLIFY=FALSE))
    
    group_candidates_list <- lapply(group_candidates_list,
                                    function(x) c(x, Res_Group =paste0("res_", x$Resolution, "_", x$group)))
    
  }
  
  
  # generate data frame
  groups_data_frame <- c()
  
  for(m in 1:M) {
    
    df_temp <- data.frame(groups[, c(1, m)])
    colnames(df_temp) <- c("SNP", "group_number")
    df_temp$Resolution <- m
    df_temp$group <- paste0("group_", df_temp$group_number)
    df_temp$res_Group <- paste0("res_",df_temp$Resolution, "_", df_temp$group)
    
    groups_data_frame <- rbind(groups_data_frame, df_temp)
    
  }
  
  return(list(groups, group_candidates_list, groups_data_frame))
  
  
}









############## GENERATE X, Y, SIGMA  ############################


# function to generate X and Y based on beta and sigma (functions for beta and sigma below)

# INPUT: 
# n, p: number of samples and predictors 
# Sigma: covariance matrix 
# beta: matrix of beta (with 0 and 1 as entries)
# Q: number of outcoems 
# myseed: seed to be used
# use_amp: whether to use signal amplitude to simulate data 
# use_snr: whether to use SNR to simulate data 
# use_distribution_beta: whether to specify a normal distirbution with specific mean and sd for the magnitude of nonzero beta 
# individual_layer_boost: whether to add some extra signal to some non-zero betas if using use_distribution_beta
# n_individual_layer_boost: specifies the number of non-zero betas that are getting the extra signal 
# mean_beta: mean of normal distribution if using use_distribustion_beta
# sd_beta: sd of normal distribution if using use_distribution_beta 
# amp: amp to be used if using use_amp 
# SNR: SNR to be used if using use_snr
# OUTPUT: list containing
# y_sampled: Y-matrix
# X: X matrix


# uses existing beta
generate_data <- function(Sigma_Dirichlet = TRUE, 
                          block_sigma = TRUE,
                          distribution_beta = TRUE, 
                          binary_y = FALSE,
                          binary_calc_method = "default",
                          beta,
                          n,
                          p,
                          Q, 
                          k, 
                          alpha0_dirichlet = 0.2, 
                          alphasum_dirichlet = 1, 
                          max_corr = 0.99, 
                          rho = 0.5,
                          binary_y_intercept = -log(99),
                          amp,
                          mean_dist_beta = 0,
                          sd_dist_beta = 0.3,
                          myseed, 
                          groups_combined, 
                          group_size_with_individual, 
                          block_size_vector, 
                          proportion_of_ones_desired) {
  
  
  
  set.seed(myseed) 
  
  # 1. Check how to generate Sigma
  if(!Sigma_Dirichlet) {
    
    
    if(!block_sigma) {
      Sigma <- generate_Sigma(method = "ar",
                              p = p,
                              rho_values = rho)
    }
    

    if(block_sigma) {
      #block_size_vector = rep(group_size_for_block_sigma, p / group_size_for_block_sigma) # just a vector containing the block-sizes
      Sigma <- make_cov(list(cov_type = "block_AR", block_sizes = block_size_vector, rho_values = rep(rho, length(block_size_vector))))

    }
    
    etas = NULL
    X <- matrix(rnorm(n * p),n) %*% chol(Sigma)
  }
  

  
  
  if(Sigma_Dirichlet) {
    
    # Create AR(k) process, see also Spector, Janson (2022) Controlled Discovery and Localization of Signals via Bayesian Linear Programming; https://github.com/amspector100/blip_sims/blob/0bd05e20011fa11b5ab04cd865662c7382a9d8c3/blip_sims/gen_data.py#L51
    alphas_dirichlet <- rep(0, k + 1)
    alphas_dirichlet[1] <- alpha0_dirichlet
    alphas_dirichlet[2:length(alphas_dirichlet)] <- (alphasum_dirichlet - alpha0_dirichlet) / k
    rhos = rdirichlet(p-1, alphas_dirichlet)
    rhos[, 1] = max(rhos[, 1], 1 - max_corr)
    rhos = rhos / rowSums(rhos)
    rhos = sqrt(rhos)
    
    etas = matrix(0, nrow = p, ncol = p)
    etas[1, 1] = 1
    
    for(j in 2:p) {
      
      rhoend = min(j, k + 1)
      
      rev_sel_rhos = rev(rhos[j-1, 2:rhoend])
      
      for (i in seq_along(rev_sel_rhos)){
        
        etas[,j] = etas[,j] + (etas[,j - i] * rev_sel_rhos[i])
        
      }
      
      scale = sqrt((1 - rhos[j - 1, 1]^2) / sum(etas[,j]^2))
      
      etas[,j] = etas[,j]*scale
      
      # Add extra noise
      etas[j, j] = rhos[j-1, 1]
      
      
    } 
    X <- matrix(rnorm(n * p),n) %*% etas
    Sigma = etas %*% t(etas)

  }
   
  
  
  # 4. Generate Y
  
  Y = matrix(NA, nrow = n, ncol = Q)
  beta_true <- matrix(NA, nrow = p, ncol = Q)
  
  for(q in 1:Q) {
    
    # Generate Beta Magnitude: Distribution or Amplitude
    
    if(distribution_beta) {
      normal_draw <- rnorm(p, mean = mean_dist_beta, sd = sd_dist_beta)
      normal_draw[abs(normal_draw < 0.1*sd_dist_beta)] <- 0.1*sd_dist_beta*sign(normal_draw[abs(normal_draw < 0.1*sd_dist_beta)])
      beta_true[, q] <- beta[, q] * normal_draw
    } else {
      beta_true[, q]  <- amp * beta[, q] / sqrt(n)
    }

    # Check whether Y is binary or not

    if(binary_y) {
      
      if(binary_calc_method == "intercept") {
        X_intercept <- cbind(rep(1, n), X)
        beta_true_intercept <- c(binary_y_intercept, beta_true[, q])
        y.sample <- function(x) rbinom(1, 1, exp(x %*% beta_true_intercept) / (1 + exp(x %*% beta_true_intercept)))
        Y[, q] <- apply(X_intercept, 1, y.sample)
      }
      
      
      if(binary_calc_method == "default"){
        beta_true[, q] <- - abs(beta_true[, q])
        y.sample <- function(x) rbinom(1, 1, ((exp(x %*% beta_true[, q]) / (1 + exp(x %*% beta_true[, q]))) - 0.5))
        Y[, q] <- apply(X, 1, y.sample)
      }
      
      if(binary_calc_method == "sparse") {
        
        prop_of_ones <- qnorm(proportion_of_ones_desired, 0, sqrt(t(beta_true[, q]) %*% Sigma %*% beta_true[, q] + pi ^ 2 / 3))
        Y[, q] <- rbinom(n, 1, plogis(prop_of_ones + X %*% beta_true[, q]))
        
      }
      
      
    }  else {
      Y[, q] <- X %*% beta_true + rnorm(n)
    }
    
   
  }
  
 
  
  if(!binary_y) {
    my_snr <- var(X %*% beta_true) 
    my_snr_alt <- var(X %*% beta_true) / (1 + var(X %*% beta_true))
  }
  
  if(binary_y) {
    my_snr <- NULL
    my_snr_alt <- NULL
  }
  
  
  
  return(list(Y = Y, X = X, etas = etas, Sigma = Sigma, beta_constructed_true = beta_true, my_snr = my_snr, my_snr_alt = my_snr_alt))
  
  
  
}



############## GENERATE BETA ############################

# function to generate beta matrix for multiple outcomes where the structure of the betas determines which outcomes are related

# INPUT: 
# p: number of predictors 
# Q: Number of outcomes
# M: number of resolutions 
# groups: data frame of dimension p x M containing the groups for each i in p for every M 
# number_non_null_big_groups: number of groups with the largest size that are non-null 
# overlap_pct: What pct of the SNPs that are non-zero should affect all outcomes? 
# nonzero_per_outcome: How many SNPs should be non-zero per outcome? 
# myseed: seed to be used 
# OUTPUT: 
# beta: matrix of dimension p x Q (with entries either 0 or 1, where a 1 indicates a non-zero effect)

# OTHER COMMENTS:
# makes sure that we sample at least one SNP for each non-zero group
# assigns equal number of non-zero snps for each group, randomly choose which groups get extra


# THIS IS USED FOR EMLKF AND UKB SIMULATION!!!! 
generate_beta_multiple_outcomes <- function(p, Q, M, groups, number_non_null_big_groups,
                                            overlap_pct, nonzero_per_outcome, myseed = 2022) {
  
  set.seed(myseed)
  
  sampled <- c()
  
  # generate matrix holding non-zero SNPs
  nonzero_matrix <- matrix(NA, nrow = nonzero_per_outcome, ncol = Q)
  
  
  # considered predictors 
  # can consider all p
  if(M == 1) {
    
    p_seq <- 1:p
    sampled_non_null_big_groups <- 1:p # can pick from all of them!
    
    
    # if M > 1, pick non-zero groups, then trickle down 
    # can only pick from the indices in the groups that are nonzero
  } else {
    
    sampled_non_null_big_groups <- sort(sample(unique(groups[, M]), number_non_null_big_groups))
    
    # now out of these non-null groups, select indices 
    indices_in_non_null_groups <- groups[groups[, M] %in%  sampled_non_null_big_groups, 1]
    
    p_seq <- indices_in_non_null_groups
    
  }
  
  
  if(overlap_pct > 0) {
    
    number_common <- floor(overlap_pct*nonzero_per_outcome)
    
    # need to make sure that there is at least one non-zero per non-zero block 
    if(nonzero_per_outcome < number_non_null_big_groups) {
      stop("Number of non-zero SNPs must be at least as large as the number of non-zero groups.")
    }
    
    # sample at least one random element per group
    # 1. equally distribute 
    # 2. randomly choose a group with gets the remainder 
    
    # how many non-zero elements per group (and how many cannot be distributed equally across groups)
    number_nonzero_per_group <- floor(number_common / number_non_null_big_groups)
    remaining_number_nonzero <- number_common %% number_non_null_big_groups
    # randomly choose groups which gets the remainder
    groups_getting_remainder <- sample(sampled_non_null_big_groups, remaining_number_nonzero)
    
    
    
    
    # sample number_nonzero_per_group randomly for each of the randomly sampled non-zero groups 
    # then randomly sample one non-zero element each out of those indices that have not been chosen previously and 
    # only for those groups that have been selected to get the remainder
    groups_df <- data.frame(groups)
    colnames(groups_df) <- paste0("layer_", seq(1, M))
    
    randomly_sampled_element_per_group <- groups_df %>%
      filter(.data[[paste0("layer_", M)]] %in% sampled_non_null_big_groups) %>%
      group_by(.data[[paste0("layer_", M)]]) %>% 
      sample_n(number_nonzero_per_group)
    
    if(length(groups_getting_remainder) > 0) {
      remainder_indices <- groups_df %>%
        filter(.data[[paste0("layer_", M)]] %in% groups_getting_remainder, 
               !(layer_1 %in% randomly_sampled_element_per_group$layer_1)) %>% 
        group_by(.data[[paste0("layer_", M)]]) %>%
        sample_n(1)
      
      
      common_nonzero <- c(remainder_indices$layer_1, randomly_sampled_element_per_group$layer_1)
      
    } else {
      common_nonzero <- c(randomly_sampled_element_per_group$layer_1)
    }
    
    
  }
  
  
  # generate non-zero for each outcome separately
  # for M > 1, this is only considering the "available" indices (i.e. those belonging to the non-null big groups)
  for(q in 1:Q) {
    
    if(overlap_pct > 0) {
      
      # if overlap_pct > 0, give common non-zero SNPs to outcome q
      nonzero_matrix[1:number_common, q] <- common_nonzero
      
      # if still non-zero SNPs that are individual for the particular outcome remain, sample them
      if(overlap_pct < 1) {
        
        available <- p_seq[!(p_seq %in% c(common_nonzero, sampled))]
        
        sampled_indiv_q <- sample(p_seq[(p_seq %in% available)], 
                                  (nonzero_per_outcome - number_common), replace = FALSE)
        
        nonzero_matrix[(number_common+1):nonzero_per_outcome, q] <- sampled_indiv_q
        sampled <- c(sampled, sampled_indiv_q)
      }
      
    }  else {
      available <- p_seq[!(p_seq %in% sampled)]
      nonzero_matrix[, q] <- sample(available, nonzero_per_outcome, replace = FALSE)
      sampled <- c(sampled, nonzero_matrix[, q])
    }
    
    
  }
  
  # put nonzero elements into matrix
  beta_ind_nonzero <- matrix(0, nrow = p, ncol = Q)
  
  for(q in 1:Q) {
    beta_ind_nonzero[groups[, 1] %in% nonzero_matrix[, q], q] <- 1
  }
  
  # randomly decide which ones are positive and which are negative 
  # this says that a single SNP is either only positive or negative for any outcome
  beta_ind_random_sign <- beta_ind_nonzero * sign(rnorm(p))
  
  return(list(beta_ind_nonzero = beta_ind_nonzero, 
              beta_ind_random_sign = beta_ind_random_sign, 
              p_seq = p_seq,
              sampled = sampled,
              sampled_non_null_big_groups = sampled_non_null_big_groups, 
              nonzero_matrix = nonzero_matrix))
  
}


generate_beta_hierarchical_with_leaf_overlap <- function(p, 
                                                         sparsity, 
                                                         overlap_pct = c(0.25, 0.25), 
                                                         Q, 
                                                         overlap_groups = list(c(1, 2), c(3, 4)), 
                                                         myseed) {
  
  
  set.seed(myseed) 
  
  number_nonzero_each = round(sparsity*p)
  
  available <- seq(1, p)
  
  
  nonzero_matrix <- matrix(NA, nrow = number_nonzero_each, ncol = Q)
  
  sampled <- c()
  counter = 1
  for(g in overlap_groups) {
    
    number_common = round(overlap_pct[counter] * number_nonzero_each)
    remaining_nonzero = number_nonzero_each - number_common
    
    # sample from available snps
    sampled_common <- sample(available, number_common)
    
    # update available snps
    sampled <- c(sampled, sampled_common)
    available <- available[!(available %in% sampled)]
    
    nonzero_matrix[1:number_common,g] = sampled_common
    
    for(i in g) {
      
      # sample remaining
      sampled_unique <- sample(available, remaining_nonzero)
      
      # update available snps
      sampled <- c(sampled, sampled_unique)
      available <- available[!(available %in% sampled)]
      
      nonzero_matrix[(number_common + 1):number_nonzero_each,i] = sampled_unique
      
    }
    
    
    
  }
  
  
  # put nonzero elements into matrix
  beta_ind_nonzero <- matrix(0, nrow = p, ncol = Q)
  
  for(q in 1:Q) {
    beta_ind_nonzero[nonzero_matrix[, q], q] <- 1
  }
  
  
  return(list(beta_ind_nonzero = beta_ind_nonzero, nonzero_matrix = nonzero_matrix))
  
  
  
  
}



random_beta_generation <- function(p, Q, sparsity, overlap_pct, myseed) {
  
  set.seed(myseed)
  
  # define number of common / individual nonzero features for each outcome
  available_features <- seq(1, p)
  
  number_nonzero = round(sparsity*p)
  number_common_nonzero = round(overlap_pct*number_nonzero)
  number_individual_nonzero = number_nonzero - number_common_nonzero
  
  if(p - number_common_nonzero - Q*number_individual_nonzero < 0) {
    stop("Increase overlap or increase sparsity. ")
  }
  
  # sample common features
  if(number_common_nonzero > 0) {
    
    sampled_common <- sample(available_features, number_common_nonzero)
    available_features <- available_features[!(available_features %in% sampled_common)]
    
  } 
  
  # generate matrix holding non-zero SNPs
  nonzero_matrix <- matrix(NA, nrow = number_nonzero, ncol = Q)
  
  # assign common nonzero features and sample individually nonzero features
  for(q in 1:Q) {
    
    if(number_common_nonzero > 0) {
      nonzero_matrix[1:number_common_nonzero, q] <- sampled_common
    }
    
    if(number_individual_nonzero > 0) {
      sampled_individual_q <- sample(available_features, number_individual_nonzero)
      available_features <- available_features[!(available_features %in% sampled_individual_q)]
      nonzero_matrix[(number_common_nonzero + 1):number_nonzero, q] <- sampled_individual_q
    }
    
    
    
    
  }
  
  
  # create indicator matrix for being nonzero based on sampled features
  beta_ind_nonzero <- matrix(0, nrow = p, ncol = Q)
  nonzero_per_q <- list()
  
  for(q in 1:Q) {
    beta_ind_nonzero[nonzero_matrix[, q], q] <- 1
    nonzero_per_q[[q]] <- nonzero_matrix[, q]
  }
  
  # randomly decide which ones are positive and which are negative 
  # this says that a single SNP is either only positive or negative for any outcome
  beta_ind_random_sign <- beta_ind_nonzero * sign(rnorm(p))
  
  nonzero_any_indiv <- unique(unlist(nonzero_per_q))
  
  return(list(nonzero_matrix = nonzero_matrix, 
              beta_ind_random_sign = beta_ind_random_sign,
              beta_ind_nonzero = beta_ind_nonzero, 
              nonzero_per_q = nonzero_per_q, 
              nonzero_any_indiv = nonzero_any_indiv))
  
  
}


get_nonzero_groups_per_q <- function(nonzero_per_q, group_df) {
  
  
  nonzero_groups_per_q <- list()
  
  for(q in 1:length(nonzero_per_q)) {
    
    nonzero_groups_per_q[[q]] <- group_df %>% dplyr::filter(individual %in% nonzero_per_q[[q]]) %>% pull(group) %>% unique()
    
    
  }
  
  return(nonzero_groups_per_q)
}


############## GENERATE BLOCK SIGMA ##############################

# Code adapted from Katsevich, Sabatti, Bogomolov (2021) Filtering the rejection set while preserving false discovery rate control

# use covariance design methods from 
# https://github.com/ekatsevi/Focused-BH/blob/c78954d3cf081fbcf22d3d7c2d26c64c17f7c7bb/src/archive/aux_make_design.R

# function to generate covariance matrix
# uses the make_cov function from https://github.com/ekatsevi/Focused-BH/blob/c78954d3cf081fbcf22d3d7c2d26c64c17f7c7bb/src/archive/aux_make_design.R

# INPUT: 
# method: whether to generate covariance matrix with block structure or using AR process 
# p: number of predictors 
# rho_values: vector containing two elements: a smaller rho, and a larger rho. 
# pct_block_rho of the nonzero groups will get the larger rho value, the rest will get the smaller rho value
# block_sizes: size of block for which block-sigma should be applied (e.g. 10)
# group_size_with_individual: vector of group sizes, e.g. c(1, 5, 10) if group sizes of 1, 5, and 10 
# nonzero_groups_per_m: list containing non-zero indices for each layer
# OUTPUT: 
# Sigma: Covariance matrix
generate_Sigma <- function(method, p, rho_values, pct_block_rho, 
                           block_sizes, group_size_with_individual, nonzero_groups_per_m) {
  
  # choose pct_block_rho of non-zero SNPs to have large rho, where block size is determined by block_sizes
  if(method == "block") {
    
    # gets indices of non-zero groups for largest group
    nonzero_block_level = nonzero_groups_per_m[[which(group_size_with_individual == block_sizes)]] 
    
    # sample pct_block_rho of those non-zero groups
    sampled_larger_rho_groups <- sample(nonzero_block_level, round(pct_block_rho*length(nonzero_block_level)))
    
    # generate vector of rho-values for each block
    # the smallest rho-value gets assigned to all blocks, except pct_block_rho of the non-zero ones
    rho_values_vector <- rep(min(rho_values), p / block_sizes)
    rho_values_vector[sampled_larger_rho_groups] <- max(rho_values)
    
    #  construct sigma 
    block_size_vector = rep(block_sizes, p / block_sizes) # just a vector containing the block-sizes
    Sigma <- make_cov(list(cov_type = "block_AR", block_sizes = block_size_vector, rho_values = rho_values_vector))
    
    
  } 
  
  # classic AR process
  if(method == "ar"){
    Sigma <- toeplitz(as.numeric(rho_values)^(0:(p-1)))
  }
  
  return(Sigma)
  
}

# See also https://github.com/ekatsevi/Focused-BH/blob/c78954d3cf081fbcf22d3d7c2d26c64c17f7c7bb/src/archive/aux_make_design.R

make_cov = function(design_options){
  stopifnot(!is.null(design_options$cov_type))
  cov_type = design_options$cov_type
  stopifnot(cov_type %in% c("AR", "block_AR"))
  if(cov_type == "AR"){
    stopifnot(!is.null(design_options$rho))
    rho = design_options$rho
    Sigma = make_AR_cov(rho, m)
  }
  if(cov_type == "block_AR"){
    block_sizes = design_options$block_sizes
    rho_values = design_options$rho_values
    Sigma = make_block_AR_cov(block_sizes, rho_values)
  }
  return(Sigma)
}

make_AR_cov = function(rho, n){
  Sigma = rho^(abs(outer(1:n, 1:n, "-")))
  return(Sigma)
}

make_block_AR_cov = function(block_sizes, rho_values){
  n = sum(block_sizes)  
  num_blocks = length(block_sizes)
  Sigma = matrix(0, n, n)
  idx = 1
  for(block in 1:num_blocks){
    block_size = block_sizes[block]
    rho = rho_values[block]
    block_idx = idx:(idx +block_size-1)
    Sigma[block_idx, block_idx] = make_AR_cov(rho, block_size)
    idx = idx + block_size
  }
  return(Sigma)
}


