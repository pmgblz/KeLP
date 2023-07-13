

              # READ AND PROCESS UKB RESULTS #

# IMPORTANT: The UKB data is not available; we have applied for it from a source (the UK Biobank) and interested parties can also apply to the same source).

# This script contains a function that reads in the lasso importance statistics
# and group information, calculates e-values, knockoff rejections / global null rejections
# for a single level of resolution


# INPUT: 
    #resolution: resolution of interest (integer)
    #population: UKB population to be studied (character)
    #alpha: alpha used for FDR control (numeric)
    #gamma: gamma used for e-value construction (numeric)
    #mydir: main directory (character) 
    #dirpath_group: directory with group participation information (character)
    #dirpath_lasso: directory with lasso results (character) 
    #outdir_path = directory where output should be saved (character) 
    #outcomes: vector containing outcome names to be considered
    #chr: chromosome of interest (integer)
    #smallest_possible_gamma: indicator for whether the smallest possible gamma should be used (leading to at least one rejection)
    #gamma_resolution_adaptive: indicator for whether a different gamma is used in different levels of resolution, 
    #use_multiplier: indicator for whether a different multiplier than the usual e-value definition should be used, 
    #multiplier = vector containing multiplier, one for each resolution, if using multiplier 
    #match_matteo_group: whether group counting should match the counting in the publicly posted results https://msesia.github.io/knockoffgwas/ukbiobank.html  
    #u: partial conjunction u (u = 1 for global)
# OUTPUT 
  # cand_groups_resolution: list, where each list element is a list for a particular group containing group-specific information (e-value, fraction, etc.) 
  # snp_group_chr_matrix_evals: matrix containing the information in cand_groups_resolution on snp-chr-group level
  # snps: vector containing all snps
  # rejected_matrix_ebh_global_res: global null rejections for this particular level of resolution (if applicable)
  # rejected_knockoff: knockoff rejections (if applicable)
  # min_non_inf_gamma = min_non_inf_gamma
 
read_process_UKB_lasso <- function(resolution, 
                                  population,
                                  alpha, 
                                  gamma,
                                  mydir, 
                                  dirpath_group = "",
                                  dirpath_lasso = "",
                                  outdir_path = paste0(my_out_dir, ""), 
                                  outcomes, 
                                  chr, 
                                  smallest_possible_gamma = FALSE, 
                                  gamma_resolution_adaptive = FALSE, 
                                  use_multiplier = FALSE, 
                                  multiplier = NA, 
                                  match_matteo_group = TRUE, 
                                  u = 1) {
  
  
  
  print(paste0("RESOLUTION IS ", resolution))
  
  
  ########### STEP 1: CREATE MATRIX CONTAINING SNP-CHR-GROUP and CHR-GROUP ###############

  ###########....... PART A: SNP-GROUP-CHR LEVEL ######################

  snp_group_chr_matrix <- c()
  
  for(c in chr) {
    
    myfilename <- paste0(mydir, dirpath_group, c, "_ibd1_res", resolution, "_grp",".txt")
    myfile <- read.csv(myfilename, sep="") 
    
    myfile$CHR <- c
    
    snp_group_chr_matrix <- rbind(snp_group_chr_matrix, myfile)
    
  }
  
  # match matteo's reported groups 
  if(match_matteo_group) {
    snp_group_chr_matrix$Group <- snp_group_chr_matrix$Group + 1
  }
  
  snp_group_chr_matrix$Resolution <- resolution
  
  
  ####################....... PART B: CHR-GROUP INFORMATION ##########################

  
  # get full list of group and chromosomes for resolution
  # creates a matrix containing two columns: CHR and Group
  # lists for each chromosome all available groups
  # if resolution is not zero, this is not the same as above
  chr_group_matrix <- c()
  
  for(c in chr) {
    
    myfilename <- paste0(mydir, dirpath_group, c, "_ibd1_res", resolution, "_grp",".txt")
    myfile <- read.csv(myfilename, sep="")
    
    chr_matrix <- matrix(NA, nrow = max(myfile$Group) + 1, ncol = 2)
    chr_matrix[, 1] <- c
    chr_matrix[, 2] <- seq(min(myfile$Group), max(myfile$Group), by = 1)
    
    chr_group_matrix <- rbind(chr_group_matrix, chr_matrix)
    
  }
  
  chr_group_matrix <- data.frame(chr_group_matrix)
  colnames(chr_group_matrix) <- c("CHR", "Group")
  
  # match matteo's reported results
  if(match_matteo_group) {
    chr_group_matrix <- chr_group_matrix %>% mutate(Group = Group + 1)
  }
  
  
  # calculate total number of groups
  total_number_groups <- dim(chr_group_matrix)[1]
  
  
  ########################## STEP 2: CALCULATE E-VALUES #######################

  # E-values are created based on chr-group level 
  # however, we do need snp-chr-group level information to identify which snps are in which groups as well!
  
  taus <- matrix(NA, nrow = 1, ncol = length(outcomes)) # taus for knockoff
  taus_eval <- matrix(NA, nrow = 1, ncol = length(outcomes)) # taus for e-values
  min_non_inf_gamma <- matrix(NA, nrow = 1, ncol = length(outcomes)) # taus for knockoff
  
  E <- matrix(NA, nrow = total_number_groups, ncol = length(outcomes) + 5)
  E[, 1] <- chr_group_matrix$CHR
  E[, 2] <- chr_group_matrix$Group
  E[, 3] <- rep(resolution, total_number_groups)
  
  
  
  fracs <- matrix(NA, nrow = total_number_groups, ncol = length(outcomes) + 5)
  fracs[, 1] <- chr_group_matrix$CHR
  fracs[, 2] <- chr_group_matrix$Group
  fracs[, 3] <- rep(resolution, total_number_groups)
  
  
  rejected_knockoff <- matrix(NA, nrow = total_number_groups, ncol = length(outcomes) + 5)
  rejected_knockoff[, 1] <- chr_group_matrix$CHR
  rejected_knockoff[, 2] <- chr_group_matrix$Group
  rejected_knockoff[, 3] <- rep(resolution, total_number_groups)
  
  # CALCULATE E-VALUE FOR EACH OUTCOME
  for(i in 1:length(outcomes)) {
    
    # files contain chr, group, W (for those with nonzero W)
    myfilename <- paste0(mydir, dirpath_lasso,"lasso_", outcomes[i], "_", population, "_res", resolution, ".txt")
    
    myfile <- read.csv(myfilename, sep="")  %>%
      dplyr::select(CHR, Group, W, BP.min, BP.max) %>%
      arrange(CHR, Group)
    
    
    # original result files only contain non-zero W. However, to calculate global null e-values, need to include those that are 
    # actually zero, so need to merge with the file containing all chr and group information
    chr_group_w_per_outcome <- dplyr::left_join(chr_group_matrix, myfile, by = c("CHR", "Group")) %>% arrange(CHR, Group)
    
    
    
    if(!gamma_resolution_adaptive) {
      
      if(smallest_possible_gamma | gamma == 0) {
        taus_eval[, i] <-  min.alpha.knockoff.threshold(myfile$W, offset = 1)[[1]]
      } else {

        
        if(resolution >= 2) {
          taus_eval[, i] <-  knockoff.threshold(myfile$W, fdr = gamma, offset = 1)
        }
        
        if(resolution == 0) {
            taus_eval[, i]  <- knockoff.threshold(myfile$W, fdr = alpha, offset = 1)
        }
        
        if(resolution == 1) {
          taus_eval[, i]  <- knockoff.threshold(myfile$W, fdr = alpha / 2, offset = 1)
        }
        
        print(paste0("GAMMA IS ", gamma))
        print(paste0("TAU IS" ,  taus_eval[, i]))
        
        
      }
      
    } 
    
    
    if(gamma_resolution_adaptive)  {

      
      if(resolution >= 2) {
        taus_eval[, i] <-  knockoff.threshold(myfile$W, fdr = gamma[resolution + 1], offset = 1)
      }
      
      if(resolution == 0) {
        
        taus_eval[, i]  <- knockoff.threshold(myfile$W, fdr = alpha, offset = 1)
        
      }
      
      
      if(resolution == 1) {
        taus_eval[, i]  <- knockoff.threshold(myfile$W, fdr = alpha / 2, offset = 1)
      }
      
      
    }
    
    min_non_inf_gamma[, i] <- min.alpha.knockoff.threshold(myfile$W, offset = 1)[[2]]
    
    print("THE SMALLEST POSSIBLE ALPHA WOULD HAVE BEEN: ")
    print(min_non_inf_gamma[, i])
    
    print(taus_eval[, i])
    
    # calculate tau for knockoff (using alpha)
    taus[, i] <-  knockoff.threshold(myfile$W, fdr = alpha, offset = 1)
    
    
    chr_group_w_per_outcome <- chr_group_w_per_outcome %>%
      mutate(W = ifelse(is.na(W), 0, W)) %>% # missing W are zero actually
      mutate(rejected_by_knockoff := ifelse(W > taus[, i], 1, 0))  # to match matteo
    
    evalue_denominator <- 1 + sum(chr_group_w_per_outcome$W <= -taus_eval[, i])
    
    if(!use_multiplier) {
      chr_group_w_per_outcome <- chr_group_w_per_outcome %>% 
        mutate(E = total_number_groups * ( (W >= taus_eval[, i]) / evalue_denominator), 
               frac = ( (W >= taus_eval[, i]) / evalue_denominator), 
               resolution = resolution) 
    } else {
      chr_group_w_per_outcome <- chr_group_w_per_outcome %>% 
        mutate(E_orig = total_number_groups * ( (W >= taus_eval[, i]) / evalue_denominator), 
               E = multiplier[resolution + 1] * ( (W >= taus_eval[, i]) / evalue_denominator),
               frac = ( (W >= taus_eval[, i]) / evalue_denominator), 
               resolution = resolution) 
    }
    
    chr_group_w_per_outcome <- chr_group_w_per_outcome %>% 
      mutate(E = total_number_groups * ( (W >= taus_eval[, i]) / evalue_denominator), 
             frac = ( (W >= taus_eval[, i]) / evalue_denominator), 
             resolution = resolution) 
    
    E[, 4] <- ifelse(is.na(E[, 4]), chr_group_w_per_outcome$BP.min, E[, 4])
    E[, 5] <- ifelse(is.na(E[, 5]), chr_group_w_per_outcome$BP.max, E[, 5])
    E[, 5+i] <- chr_group_w_per_outcome$E
    
    fracs[, 4] <- ifelse(is.na(fracs[, 4]), chr_group_w_per_outcome$BP.min, fracs[, 4])
    fracs[, 5] <- ifelse(is.na(fracs[, 5]), chr_group_w_per_outcome$BP.max, fracs[, 5])
    fracs[, 5+i] <- chr_group_w_per_outcome$frac
    
    rejected_knockoff[, 4] <- ifelse(is.na(rejected_knockoff[, 5]), chr_group_w_per_outcome$BP.min, rejected_knockoff[, 5])
    rejected_knockoff[, 5] <- ifelse(is.na(rejected_knockoff[, 5]), chr_group_w_per_outcome$BP.max, rejected_knockoff[, 5])
    rejected_knockoff[, 5+i] <- chr_group_w_per_outcome$rejected_by_knockoff
    
    
    # save file
    write.csv(chr_group_w_per_outcome, paste0(outdir_path, "snp_group_w_",
                                              outcomes[i], "_" ,population ,"_res", resolution, ".txt"), row.names = FALSE)
    
    
  }
  
  
  # E VALUES
  E_means_matrix <- matrix(NA, nrow = total_number_groups, ncol = 6)
  E_means_matrix[, 1:5] <- E[, 1:5]
  
  # calculate means
  if(length(outcomes) == 1) {
    E_means_matrix[, 6] <- E[, 6]
  } else {
    E_means_matrix[, 6]  <- partial_conjunction_averaging(E[, 6:(dim(E)[2])], Q = length(outcomes), u = u)$avg
  }
  
  
  colnames(E_means_matrix) <- c("CHR", "Group","Resolution","BP.min" , "BP.max","evals")
  E_means_matrix <- as.data.frame(E_means_matrix)
  
  
  # FRACTIONS
  fracs_means_matrix <- matrix(NA, nrow = total_number_groups, ncol = 6)
  fracs_means_matrix[, 1:5] <- E[, 1:5]
  
  # calculate means
  if(length(outcomes) == 1) {
    fracs_means_matrix[, 6] <- fracs[, 6]
  } else {
    fracs_means_matrix[, 6] <- partial_conjunction_averaging(fracs[, 6:(dim(fracs)[2])], Q = length(outcomes), u = u)$avg
  }
  
  
  
  colnames(fracs_means_matrix) <- c("CHR", "Group","Resolution","BP.min" , "BP.max","fracs")
  fracs_means_matrix <- as.data.frame(fracs_means_matrix)
  
  
  snp_group_chr_matrix_evals <- dplyr::left_join(snp_group_chr_matrix, E_means_matrix, by = c("CHR", "Group", "Resolution"))
  snp_group_chr_matrix_evals <- dplyr::left_join(snp_group_chr_matrix_evals, 
                                                 fracs_means_matrix %>% dplyr::select(Group, CHR, Resolution, fracs), by = c("CHR", "Group", "Resolution"))
  print(sum(is.na(snp_group_chr_matrix_evals$evals)))
  print(dim(snp_group_chr_matrix))
  print(dim(snp_group_chr_matrix_evals))
  
  # Now create list
  snp_group_chr_matrix_evals$CHR_Group <- paste0(snp_group_chr_matrix_evals$CHR, "_", snp_group_chr_matrix_evals$Group)
  
  
  # it seems that even in higher resolutions there sometimes are groups with only one SNP
  cand_groups_resolution <- split(snp_group_chr_matrix_evals, snp_group_chr_matrix_evals$CHR_Group)
  cand_groups_resolution <- lapply(cand_groups_resolution, function(x) create_list_element(list_element = x))
  cand_groups_resolution <- unname(cand_groups_resolution)
  
  
  
  # TEST GLOBAL NULL BASED ON EVALUES
  E_means <- fracs_means_matrix[, 6]*total_number_groups # for rejecting global null need to use normal e-values to run ebh!
  E_means_ord <- order(E_means, decreasing = TRUE)
  E_means <- sort(E_means, decreasing = TRUE)
  comp_global <- E_means >= (total_number_groups / alpha / (1:total_number_groups))
  
  id_global <- max(which(comp_global > 0))
  if(id_global > 0){
    rej_global <- E_means_ord[1:id_global]
  }else{
    rej_global<- NULL
  }
  
  rejected_matrix_ebh_global_res <- E_means_matrix[rej_global, ]
  
  
  
  if(length(outcomes) == 1) {
    return(list(cand_groups_resolution = cand_groups_resolution, 
                snp_group_chr_matrix_evals = snp_group_chr_matrix_evals,
                snps = unique(snp_group_chr_matrix$SNP), 
                rejected_matrix_ebh_global_res = rejected_matrix_ebh_global_res, 
                rejected_knockoff = rejected_knockoff, 
                min_non_inf_gamma = min_non_inf_gamma))
  } else {
    return(list(cand_groups_resolution = cand_groups_resolution, 
                snp_group_chr_matrix_evals = snp_group_chr_matrix_evals,
                snps = unique(snp_group_chr_matrix$SNP), 
                rejected_matrix_ebh_global_res = rejected_matrix_ebh_global_res, 
                min_non_inf_gamma = min_non_inf_gamma))
  }
  
  
  
}


# helper function to create a single list element (that is a list)
# INPUT: 
  # list_element: element to be converted to list 
# OUTPUT: 
  # x_list: converted list element as list

create_list_element <- function(list_element) {
  
  x_list <- as.list(list_element)
  x_list <- lapply(x_list, function(x) unique(x))
  
  return(x_list)
  
  
}

