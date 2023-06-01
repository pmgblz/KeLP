

# e-BH related functions #


#1. function for running e-BH 
#2. functions for running focused-eBH 


########################## e-BH ######################################



# INPUT: 
# E: vector of e-values or fractions 
# alpha: alpha for fdr-control (numeric)
# number_groups: total number of groups (integer)
# use_fractions: indicator for whether to run e-BH with on fractions (e-value without the multiplier of the total number of groups) or on e-values

# OUTPUT: 
# rej: rejections

run_eBH <- function(E, alpha, number_groups, use_fractions = TRUE) {
  
  
  E_ord <- order(E, decreasing = TRUE)
  E <- sort(E, decreasing = TRUE)
  
  if(use_fractions) {
    comp <- E >= (1 / alpha / (1:number_groups)) # run on fractions (e-values not multiplier by p, so here also no p (equivalent to running on original e-values))
  } else {
    comp <- E >= (number_groups / alpha / (1:number_groups)) # run on e-values multiplied by number_groups
  }
  
  id <- max(which(comp > 0))
  
  if(id > 0){
    rej <- E_ord[1:id]
  } else{
    rej <- NA
  }
  
  return(rej)
  
}



######################## FOCUSED E-BH ####################

# Function to run focused e-BH in simulation setting. 

# INTPUT: 
# fracs: list of fractions (i.e. e-valuves without the multiplier), length of list = M (number of layers)
# M: number of layers 
# group_candidates_list: list of length = number of groups, where each list element contains group information, including the fraction, resolution, etc. 
# group_data_frame: contains group information as a data frame, created in sim_data_aux() 
# numerator_threshold: numerator of self-consistency threshold (i.e. M if using fractions in multiple layers)
# alpha: desired level of FDR control 


# OUTPUT: 
# R_final_outer_nodes: rejections by focused e-BH 

focusedeBH <- function(fracs, M, group_candidates_list, group_data_frame, numerator_threshold, alpha = 0.2) {
  
  ts_orig = unlist(fracs)
  ts_unique <- c(unique(ts_orig), Inf)
  
  trt = c()
  all_resolutions <- seq(1, M)
  
  # compute rejections for all possible thresholds
  for(i in 1:length(ts_unique)) {
    
    # get rejections 
    R = which(ts_orig >= ts_unique[i])
    R_list_elements <- group_candidates_list[R]
    
    if(length(R_list_elements) > 0) {
      
      # only keep those with maximum fraction
      R_elements <- plyr::ldply(R_list_elements, data.frame) %>% 
        group_by(Res_Group) %>% 
        dplyr::filter(fracs == max(fracs))
      
      # filter to outer nodes
      R_elements_info <- group_data_frame %>% 
        dplyr::filter(res_Group %in% unique(R_elements$Res_Group)) %>% 
        mutate(Resolution = Resolution)
      
      R_outer_nodes <- filter_outer_node(R_elements_info,
                                         resolutions = all_resolutions,
                                         group_identifiyer = "res_Group")
    } else {
      R_outer_nodes <- NULL
    }
    
    if(!is.null(R_outer_nodes)) {
      trt[i] <- ts_unique[i]*nrow(R_outer_nodes)
    } else {
      trt[i] <- 0
    }
    
    
  }
  
  # compute final threshold
  final_treshold <- min(ts_unique[which(trt >= numerator_threshold / alpha)])
  
  # get final rejections 
  R_final = which(ts_orig >= final_treshold)
  R_final_list_elements <- group_candidates_list[R_final]
  
  if(length(R_final_list_elements) > 0) {
    
    # only keep those with maximum fraction
    R_final_elements <- plyr::ldply(R_final_list_elements, data.frame) %>% 
      group_by(Res_Group) %>% 
      dplyr::filter(fracs == max(fracs))
    
    # filter to outer nodes
    R_final_elements_info <- group_data_frame %>% 
      dplyr::filter(res_Group %in% unique(R_final_elements$Res_Group)) 
    
    R_final_outer_nodes <- filter_outer_node(R_final_elements_info,
                                             resolutions = all_resolutions,
                                             group_identifiyer = "res_Group")
    
  } else {
    R_final_outer_nodes <- NULL
  }
  
  
  return(R_final_outer_nodes)
  
}



# Function to run focused e-BH in UKB setting. 

# INPUT: 
# unique_groups: data frame containing group-level information and fractions / e-values
# snp_groups: data frame containing information on in which group each snp is belonging
# numerator_threshold: numerator of self-consistency threshold (i.e. M if using fractions in multiple layers)
# alpha: desired level of FDR control


# OUTPUT: 
# R_final_outer_nodes: rejections by focused e-BH 

focusedeBH_UKB <- function(unique_groups, snp_groups, numerator_threshold, alpha = 0.1) {
  
  
  set.seed(2022)
  
  ts_orig <- unique_groups$fracs
  ts_unique <- c(unique(ts_orig), Inf)
  ts_unique <- ts_unique[!(ts_unique == 0)]
  
  trt = c()
  all_resolutions <- seq(0, 6)
  
  # compute rejections for all possible thresholds
  for(i in 1:length(ts_unique)) {
    
    print(i)
    
    # get rejections 
    R = which(ts_orig >= ts_unique[i])
    
    # if duplicate groups, choose the one with highest fraction 
    R_elements = unique_groups[R, ]
    
    
    if(dim(R_elements)[1] > 0) {
      R_elements_unique <- R_elements %>% 
        group_by(Res_CHR_Group) %>% 
        dplyr::filter(fracs == max(fracs)) # groups can appear multiple times if having multiple outcomes, select the group with the higher fraction
    } else {
      R_elements_unique <- R_elements
    }
    
    # get snps corresponding to selected groups 
    snp_groups_matrix <- snp_groups %>% dplyr::filter(Res_CHR_Group %in% R_elements_unique$Res_CHR_Group)
    
    # filter to outer nodes
    R_outer_nodes <- filter_outer_node(snp_groups_matrix,
                                       resolutions = all_resolutions,
                                       group_identifiyer = "Res_CHR_Group")
    
    if(!is.null(R_outer_nodes)) {
      trt[i] <- ts_unique[i]*nrow(R_outer_nodes)
    } else {
      trt[i] <- 0
    }
    
    
  }
  
  # compute final threshold
  final_treshold <- min(ts_unique[which(trt >= numerator_threshold / alpha)]) 
  
  # get rejection set 
  R_final = which(ts_orig >= final_treshold)
  
  # if duplicate groups, choose the one with highest fraction 
  R_final_elements = unique_groups[R_final, ]
  
  
  if(dim(R_elements)[1] > 0) {
    R_final_elements_unique <- R_final_elements %>% 
      group_by(Res_CHR_Group) %>% 
      dplyr::filter(fracs == max(fracs)) # groups can appear multiple times if having multiple outcomes, select the group with the higher fraction
  } else {
    R_final_elements_unique <- R_final_elements
  }
  
  # get snps corresponding to selected groups 
  snp_groups_matrix_final <- snp_groups %>% dplyr::filter(Res_CHR_Group %in% R_final_elements_unique$Res_CHR_Group)
  
  # filter to outer nodes
  R_final_outer_nodes <- filter_outer_node(snp_groups_matrix_final,
                                           resolutions = all_resolutions,
                                           group_identifiyer = "Res_CHR_Group")
  
  
  
  R_final_outer_info <- R_final_elements_unique %>% dplyr::filter(Res_CHR_Group %in% R_final_outer_nodes$Res_CHR_Group)
  
  fracs_corresponding <- unique(R_final_outer_info$fracs)
  
  if(sum(fracs_corresponding < numerator_threshold / (alpha * nrow(R_final_outer_nodes))) > 0) {
    stop("Not self consistent!")
  }
  
  return(R_final_outer_nodes)
  
}



# Function to run focused e-BH in hierarchical outcome setting. 

# INPUT: 
# fracs: list of fractions, each sub-list is an outcome
# fracs_df: data frame containing fracions, and further info (e.g. on sub-tree, level, see hierarchical_outcomes.R for more information)
# numerator_threshold: numerator of self-consistency threshold (i.e. M if using fractions in multiple layers)
# alpha: desired level of FDR control


# OUTPUT: 
# R_final_outer_nodes: rejections by focused e-BH 


#assumes uniform multiplier
focusedeBH_hierarchical_y <- function(fracs, 
                                      fracs_df,
                                      numerator_threshold = Q,
                                      alpha = 0.2) {
  
  ts_orig = unlist(fracs)
  ts_unique <- c(unique(ts_orig), Inf)
  ts_unique <- ts_unique[!(ts_unique == 0)]
  
  trt = c()
  
  # compute rejections for all possible thresholds
  for(i in 1:length(ts_unique)) {
    
    # get rejections 
    R = which(ts_orig >= ts_unique[i])
    R_elements <- fracs_df[R, ]
    
    if(dim(R_elements)[1] > 0) {
      
      # for each "subtree" or root-node, filter to max outcome 
      # this keeps the most specific rejected node among each of the "sub-trees" (7), (1, 2, 5), (3, 4, 6)
      R_outer_nodes <- R_elements %>% 
        group_by(SNP, subtree) %>% 
        dplyr::filter(level == max(level))
      
      # for each snp, check whether another outcome was rejected that is not Y1
      max_subtree_by_snp <- R_outer_nodes %>% 
        group_by(SNP) %>% 
        mutate(max_subtree = max(subtree)) %>% 
        dplyr::select(SNP, max_subtree) %>% 
        unique()
      
      R_outer_nodes <- left_join(R_outer_nodes, max_subtree_by_snp, by = "SNP")
    
      # if Y1 was rejected but also a node further down the tree, only keep the node further down
      # however, we can keep rejected nodes in different subtrees
      R_outer_nodes <- R_outer_nodes %>% 
        filter(!(subtree == 1 & max_subtree > 1))
      
      
    } else {
      R_outer_nodes <- NULL
    }
    
    
    if(!is.null(R_outer_nodes)) {
      trt[i] <- ts_unique[i]*nrow(R_outer_nodes)
    } else {
      trt[i] <- 0
    }
    
    
  }
  
  
  # compute final threshold
  final_treshold <- min(ts_unique[which(trt >= numerator_threshold / alpha)]) 
  
  # get rejection set 
  R_final = which(ts_orig >= final_treshold)
  
  R_final_elements <- fracs_df[R_final, ]
  
  if(dim(R_final_elements)[1] > 0) {
    
    
    
    # for each "subtree" or root-node, filter to max outcome 
    # this keeps the largest rejected node among each of the "sub-trees" (1), (2, 4, 5), (3, 6, 7)
    R_final_outer_nodes <- R_final_elements %>% 
      group_by(SNP, subtree) %>% 
      dplyr::filter(level == max(level))

    
    # for each snp, check whether another outcome was rejected that is not Y1
    final_max_subtree_by_snp <- R_final_outer_nodes %>% 
      group_by(SNP) %>% 
      mutate(max_subtree = max(subtree)) %>% 
      dplyr::select(SNP, max_subtree) %>% 
      unique()
    
    R_final_outer_nodes <- left_join(R_final_outer_nodes, final_max_subtree_by_snp, by = "SNP")
    
    # if Y1 was rejected but also a node further down the tree, only keep the node further down
    # however, we can keep rejected nodes in different subtrees
    R_final_outer_nodes <- R_final_outer_nodes %>% 
      filter(!(subtree == 1 & max_subtree > 1))

    
  } else {
    R_final_outer_nodes <- NULL
  }    
  
  
  if(sum(R_final_outer_nodes$fracs < numerator_threshold / (alpha * nrow(R_final_outer_nodes))) > 0) {
    
    print(summary(R_final_outer_nodes$fracs))
    print(numerator_threshold / (alpha * nrow(R_final_outer_nodes)))
    print(sum(R_final_outer_nodes$fracs < numerator_threshold / (alpha * nrow(R_final_outer_nodes))))
    print(which(R_final_outer_nodes$fracs < numerator_threshold / (alpha * nrow(R_final_outer_nodes))))
    
    stop("Not self consistent!")
  }
  
  return(R_final_outer_nodes)
  
}

