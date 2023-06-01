
# Functions for e-MLKF #


# This file contains functions for the multilayer knockoff filter and the e-multilayer knockoff filter


############# Functions modified from Katsevich, Sabatti (2019): Multilayer Knockoff filter: Controlled variable selection at multiple resolutions and Barber, Ramdas (2017) The p-filter: multi-layer FDR control for grouped hypotheses ###############

# function to find the number of rejections in a particular element 

# INPUT: 
  # combs_element: particular combination of gammas, one for each layer
  # alpha_emlkf: alpha for FDR control in e-MLKF 
  # groups: data frame containing group structure information 
  # my_offset: offset used for knockoff 
  # G: containing group size information
  # W: importance statistics as list of length M

# OUTPUT: list containing

  # num_rejections: number of rejections 
  # taus_gamma: knockoff stopping time corresponding to the particular gamma

find_num_rejections <- function(combs_element, alpha_emlkf, groups, my_offset, G, W, myseed) {
  
  set.seed(myseed)
  
  combs_element <- as.numeric(combs_element)
  E_gamma <- list()
  taus_gamma <- c()

  
  for(m in 1:M){
    
    # no replacement of gamma
    taus_gamma[m] <- knockoff.threshold(W[[m]], fdr = as.numeric(combs_element[m]), offset = my_offset)
    

    E_gamma[[m]] <- G[m] * (W[[m]] >= taus_gamma[m]) / (1 + sum(W[[m]] <= -taus_gamma[m]))
    
  }
  
  S_hat_alpha <- efilter_group_eval(E = E_gamma, alpha_emlkf = alpha_emlkf, groups = groups)
  num_rejections <- length(S_hat_alpha)
  
  return(list(num_rejections = num_rejections, taus_gamma = taus_gamma))
}




# multilayer knockoff filter with e-values
# parts of code-structure taken from Katsevich, Sabatti (2019), downloaded from https://katsevich-lab.github.io/publications/

# INPUT: 
  # X 
  # Y
  # W: importance statistics, list of length M (number of layers)
  # groups is a n-by-M matrix; groups[i,m] = which group does E[i] belong to for the m-th grouping
  # final_gamma_combination: gammas to be used for e-value creation 
  # alpha_emlkf in  [0,1]^M = vector of target FDR levels


# OUTPUT: 

  # S_hat: rejections by layer (list of length M)

e_multilayer_knockoff_filter <- function(X, Y, W, 
                                         groups, 
                                         final_gamma_combination, 
                                         alpha_emlkf) {
  
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
  
  # compute e-values for each layer based on final gamma combination
  E <- list()
  taus <- c()
  for(m in 1:M){
    taus[m] <- knockoff.threshold(W[[m]], fdr = final_gamma_combination[m], offset = my_offset)
    E[[m]] <- G[m] * (W[[m]] >= taus[m]) / (1 + sum(W[[m]] <= -taus[m]))
  }
  
  
  # run e-filter on the e-values
  S_hat <- efilter_group_eval(E = E, alpha_emlkf = alpha_emlkf, groups = groups)
  
  return(list(S_hat = S_hat))
  
}



# e-filter
# Code modified form p-filter code in Barber, Ramdas (JRSS B), downloaded from https://rinafb.github.io/research/

# INPUT: 
  # E: e-values; list of length M
  # alpha_emlkf in  [0,1]^M = vector of target FDR levels
  # groups is a n-by-M matrix; groups[i,m] = which group does E[i] belong to for the m-th grouping

efilter_group_eval = function(E,alpha_emlkf,groups){

  n = nrow(groups)
  M = ncol(groups)
  G = apply(groups,2,max) # G[m] = # groups, for grouping m
  
  # initialize
  thresh = 1/alpha_emlkf
  Sh = 1:n
  for(m in 1:M){
    pass_means_m = which(is.element(groups[,m],which(E[[m]]>=thresh[m])))
    Sh = intersect(Sh,pass_means_m)
  }
  done = FALSE
  
  while(!done){
    thresh_old = thresh
    for(m in 1:M){
      # which groups, for the m-th grouping, 
      #	have any potential discoveries?
      Shm = sort(unique(groups[Sh,m]))
      
      # run BH, constraining to Shm
      Evals_m = rep(-0.01,G[m]); # <0 for groups not in Dm
      Evals_m[Shm] = E[[m]][Shm]
      
      khatm = max(0,which(sort(Evals_m, decreasing = TRUE)>=G[m]/((1:G[m])*alpha_emlkf[m])))
      thresh[m] = G[m]/(khatm*alpha_emlkf[m])
      
      Sh = intersect(Sh,
                     which(is.element(groups[,m],which(Evals_m>=thresh[m]))))
    }
    if(all(thresh_old==thresh)){done = TRUE}
  }
  
  return(Sh)
  
}




############### The following functions are taken directly KATSEVICH, SABATTI (2019, Annals of Applied Statistics, Multilayer Knockoff Filter) #####################

# code downloaded from: https://katsevich-lab.github.io/publications/
# See katsevich, sabatti (2019) for further details on these functions

# get estimated number of false discoveries at each layer for given thresholds
get_V_hats = function(P, thresh_idx, FDP_hat_type){
  M = length(thresh_idx)
  
  V_hats = numeric(M)
  switch(FDP_hat_type,
         "kn" = {
           V_hats[thresh_idx == 0] = 0
           for(m in which(thresh_idx > 0)){
             V_hats[m] = sum(P[[m]][1:thresh_idx[m]] == 1)
           }    
         },
         
         "kn+" = {
           V_hats[thresh_idx == 0] = 1
           for(m in which(thresh_idx > 0)){
             V_hats[m] = 1 + sum(P[[m]][1:thresh_idx[m]] == 1)
           }    
         }
  )
  return(V_hats)
}

# get selection set for given thresholds
get_S_hat = function(P, groups, thresh_idx){
  n = nrow(groups)
  M = ncol(groups)
  if(any(thresh_idx == 0)){
    S_hat = numeric(0)    
  }
  else{
    S_hat = 1:n # current selection set
    for(m in 1:M){
      S_tilde_m = which(is.element(groups[,m],intersect(which(P[[m]] == 0), 1:thresh_idx[m])))
      S_hat = intersect(S_hat, S_tilde_m)
    }
  }
  return(S_hat)
}


multilayer_knockoff_filter = function(X, Y, W, groups, alpha_mlkf, FDP_hat_type){
  
  # problem dimensions
  N = nrow(X) # number of samples
  n = nrow(groups) # number of variables
  M = ncol(groups) # number of layers
  G = apply(groups,2,max) # G[m] = number of groups for at layer m
  
  # check input for correctness
  stopifnot(length(Y) == N)
  stopifnot(ncol(X) == n)
  stopifnot(length(alpha_mlkf) == M)
  stopifnot(FDP_hat_type %in% c("kn", "kn+"))
  
  # reorder groups based on magnitudes of knockoff statistics
  group_orders = list()               # ordering of groups at each layer
  groups_reordered = matrix(0, n, M)  # reordered group assignments
  for(m in 1:M){
    group_orders[[m]] = order(abs(W[[m]]), decreasing = TRUE)
    groups_reordered[,m] = invPerm(group_orders[[m]])[groups[,m]]
  }
  
  # define one-bit p-values
  P = list()
  allowable_thresh_idx = list()
  for(m in 1:M){
    kn_stats_ordered = W[[m]][group_orders[[m]]]
    signs = sign(kn_stats_ordered)
    allowable_thresh_idx[[m]] = which(abs(kn_stats_ordered[1:(G[m]-1)]) > abs(kn_stats_ordered[2:G[m]]))
    allowable_thresh_idx[[m]] = c(allowable_thresh_idx[[m]], G[m])
    P[[m]] = rep(1, G[m])
    P[[m]][signs == 1] = 0
  }
  
  # run filter to get threshold indices at each layer
  thresh_idx = get_thresholds(P, allowable_thresh_idx, groups_reordered, alpha_mlkf, FDP_hat_type)
  
  # extract knockoff statistic thresholds from threshold indices
  thresh = numeric(M)
  for(m in 1:M){
    if(thresh_idx[m] == 0){
      thresh[m] = Inf
    }
    else{
      thresh[m] = abs(W[[m]][group_orders[[m]][thresh_idx[m]]])
    }
  }
  
  # get selection set
  S_hat = get_S_hat(P, groups_reordered, thresh_idx)
  
  # return output: selection set, threshold indices, and thresholds
  output <- list()
  output$S_hat = S_hat
  output$thresh_idx = thresh_idx
  output$thresh = thresh
  return(output)
}

# find multilayer knockoff filter thresholds at each layer
get_thresholds = function(P, allowable_thresh_idx, groups, q, FDP_hat_type){
  cat(sprintf("Searching for thresholds...\n"))
  n = nrow(groups)
  M = ncol(groups)
  
  G = apply(groups,2,max) # G[m] = # groups, for grouping m
  
  # initialize thresholds
  thresh_idx = G
  
  # find thresholds
  done = FALSE
  while(!done){
    thresh_idx_old = thresh_idx
    for(m in 1:M){
      if(thresh_idx[m] >= 1){
        thresh_idx_m = 0
        for(thresh_idx_m_tmp in thresh_idx[m]:1){
          if(thresh_idx_m_tmp %in% allowable_thresh_idx[[m]]){
            thresh_idx_tmp = thresh_idx
            thresh_idx_tmp[m] = thresh_idx_m_tmp
            FDP_hat = get_FDP_hat(P, groups, thresh_idx_tmp, FDP_hat_type)
            if(FDP_hat[m] <= q[m]){
              thresh_idx_m = thresh_idx_m_tmp
              break
            }
          }
        }
        thresh_idx[m] = thresh_idx_m
        if(thresh_idx[m] == 0){
          thresh_idx = rep(0, M)
          done = TRUE
          break
        }
      }
    }
    if(all(thresh_idx_old==thresh_idx)){done = TRUE}
  }
  
  cat("done.\n")
  return(thresh_idx)  
}


# get FDP-hat for given thresholds 
get_FDP_hat = function(P, groups, thresh_idx, FDP_hat_type){
  n = nrow(groups)
  M = ncol(groups)
  
  S_hat = get_S_hat(P, groups, thresh_idx)
  
  S_hats_m = sapply(1:M, function(m)(length(unique(groups[S_hat,m]))))    
  V_hats_m = get_V_hats(P, thresh_idx, FDP_hat_type)
  
  FDP_hat = V_hats_m/(pmax(1, S_hats_m))
  
  return(FDP_hat)
}


