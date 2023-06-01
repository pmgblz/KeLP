
# functions for running kelp, some code snippets are inspired by Spector, Janson (2022) Controlled Discovery and Localization of Signals via Bayesian Linear Programming; https://github.com/amspector100/blipr/blob/main/R/blip.R

#1. function to create location constraint matrix 
#2. function to run kelp
#3. kelp utility functions (e.g. calculating inverse size)

####################### LOCATION CONSTRAINT MATRICES ##########################


# function to create location constraint matrix
# INPUT:
  # group_candidates_list: list containing groups, each element is a list containing information about each group, see sim_data_aux
  # individual_snps_candidates: vector of individual snps 
  # divisor: (integer) number of parts (+ remainder) in which the location constraint matrix will be created  
  # n_groups_filtered: number of groups remaining after filtering out all groups with zero e-values (integer)
  # n_individual_hypotheses_filtered: number of individual hypotheses (or snps) remaining after filtering out all groups with zero e-values

# OUTPUT: 
  # A: location constraint matrix

create_location_constraint_matrix <- function(group_candidates_list,
                                              individual_snps_candidates, 
                                              divisor = 4,
                                              n_groups_filtered,
                                              n_individual_hypotheses_filtered) {
  
  
  set.seed(2022)
  
  # parts in which A will be created
  parts <- seq(1, divisor)
  
  A <- c()
  
  prev_part <- 0
  
  step_size <- floor((length(group_candidates_list) / divisor))
  remaining <- n_groups_filtered- step_size*divisor
  
  for(p in parts) {
    
    part_seq <- seq(max(prev_part) + 1, floor((length(group_candidates_list) / divisor)) * p)
    
    A_part <- matrix(0, length(part_seq), n_individual_hypotheses_filtered)
    colnames(A_part) <- paste0("SNP_", individual_snps_candidates)
    
    counter <- 0
    for (j in part_seq) {
      counter <- counter + 1
      A_part[counter, paste0("SNP_",group_candidates_list[[j]]$SNP)] = 1
    }
    
    
    A <- rbind(A, A_part)
    prev_part <- part_seq
    
  }
  
  
  # take care of remaining
  if(remaining > 0) {
    part_seq <- seq(max(prev_part) + 1, length(group_candidates_list))
    
    A_part <- matrix(0, length(part_seq), n_individual_hypotheses_filtered)
    colnames(A_part) <- paste0("SNP_", individual_snps_candidates)
    
    counter <- 0
    for (j in part_seq) {
      counter <- counter + 1
      A_part[counter, paste0("SNP_",group_candidates_list[[j]]$SNP)] = 1
    }
    
    
    A <- rbind(A, A_part)
    prev_part <- part_seq
  }
  
  
  
  return(A)
}


########################## KELP ###################


# function to run kelp, contains code snippets from https://github.com/amspector100/blipr/blob/main/R/blip.R
# INPUT:
  # group_candidates_list: list containing groups, each element is a list containing information about each group, see sim_data_aux
  # individual_snps_candidates: vector of individual snps 
  # single_u: indicator for whether a single partial conjunction u is tested 
  # alpha: alpha used for FDR control (numeric)
  # verbose: indicator for whether to print additional output / information when running optimization procedure 
  # weighted: indicator for whether using weighted optimization (weight by inverse group size)
  # fractions_only: indicator whether to run optimization based on fractions (uniform multiplier) or on original e-values 
  # M: number of resolutions (integer)
  # n_individual_hypotheses_filtered: number of individual hypotheses (or snps) remaining after filtering out all groups with zero e-values
  # n_groups_filtered: number of groups remaining after filtering out all groups with zero e-values (integer)
  # total_groups: total number of groups (integer)
  # Q: number of outcomes (integer)
  # u_seq: sequence of u's used for partial conjunction 

# OUTPUT: list containing 
  # detections: elements of group_candidates_list that were rejected
  # group_candidates_list: entire group_candidates_list, including indicator in each list-element for whether group was rejected or not

kelp <- function(group_candidates_list, 
                   individual_snps_candidates, 
                   single_u = TRUE, 
                   alpha, 
                   verbose = F,
                   weighted = T, 
                   fractions_only = TRUE,
                   M, # number of resolutions
                   n_individual_hypotheses_filtered, 
                   n_groups_filtered, 
                   total_groups, 
                   Q, # number of outcomes
                   u_seq = 1) {
  
  set.seed(2022)
  
  # this is not necessarily the same as total groups, as cand-groups can be pre-filtered
  n_groups_filtered <- length(group_candidates_list)
  
  # calculate (inverse groups size) weight for each group
  #weights <- sapply(group_candidates_list, inverse_size_alpha, alpha = alpha)
  weights <- sapply(group_candidates_list, inverse_size)
 
  # resolutions 
  inverse_resolution <- sapply(group_candidates_list, get_inverse_resolution)
  
  # get partial value
  if(!single_u) {
    u <- sapply(group_candidates_list, function(x) x$u)
  }
  
  # Extract e-values or fractions (for e-values with uniform weighting)
  if(fractions_only) {
    evals <- sapply(group_candidates_list, function(x) x$fracs)
  } else {
    evals <- sapply(group_candidates_list, function(x) x$evals)
  }
  
  
  # Constraints to ensure selected groups are disjoint
  # split into parts for faster processing
  print("Creating Location Constraint Matrix ...")
  A <- create_location_constraint_matrix(group_candidates_list, 
                                         individual_snps_candidates, 
                                         divisor = 4, 
                                         n_groups_filtered = n_groups_filtered,
                                         n_individual_hypotheses_filtered = n_individual_hypotheses_filtered)
  
  
  # Assemble variables and constraints
  # does not include e-value in objective
  print("Setting up Optimization Problem ...")
  x <- CVXR::Variable(n_groups_filtered, integer=TRUE) # already force x to be an integer, x is variable we are optimizing over
  
  if(weighted) {
    if(single_u) {
      objective <- CVXR::Maximize(sum(weights * inverse_resolution * x))
    } else {
      objective <- CVXR::Maximize(sum(weights * inverse_resolution * u * x))
    }
    
  } else {
    objective <- CVXR::Maximize(sum(x)) 
  }
  
  # constraints: only one per location
  b <- rep(1, n_individual_hypotheses_filtered)
  constraints <- list(
    x >= 0,
    x <= 1,
    t(A) %*% x <= b # makes sure that regions are disjoint
  )
  
  
  # FDR constraint
  
  # check if adjustments for alpha need to be made
  alpha_mod <- alpha / sum(Q / (Q - u_seq + 1))
  print(alpha_mod)
  
  if(fractions_only) {
    constraints <- c(constraints, list(
      M - evals*alpha_mod*sum(x) <= M*(1-x)
    ))
  } else {
    constraints <- c(constraints, list(
      total_groups - evals*alpha_mod*sum(x) <= total_groups*(1-x)
    ))
  }
  
  
  # post optimization problem
  problem <- CVXR::Problem(
    objective=objective, constraints=constraints
  )
  
  
  # solve the problem (use GLPK solver since we force x to be an integer)
  print("Solving Optimization Problem ...")
  result <- CVXR::solve(
    problem, solver='GLPK', verbose=verbose)
  
  # obtain selections
  selections <- as.numeric(result$getValue(x))
  
  # Save information
  for (j in 1:n_groups_filtered) {
    group_candidates_list[[j]]$selected <- selections[j]
  }
  
  detections <- group_candidates_list[lapply(group_candidates_list, function(x) x$selected) == 1]
  
  return(list(detections, group_candidates_list))
  
  
}


# get inverse resolution (+ 1 because resolution can start at 0)
get_inverse_resolution <- function(group_candidates_list) {
 return((1 / (group_candidates_list$Resolution + 1)))
}

# calculate inverse size of group
inverse_size <- function(group_candidates_list) {
  return(1 / (length(group_candidates_list$SNP)))
}

inverse_size_alpha <- function(group_candidates_list, alpha) {
  weight <- alpha + (1 - alpha) / (length(group_candidates_list$SNP))
}  