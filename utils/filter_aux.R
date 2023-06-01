
                  # Filtering Functions #


# This file contains functions to filter to the outer nodes and to compare two levels of resolution:


# function to compare two levels of resolution: For each (nested) group-structure spanning two levels of resolution, 
# only keep the group with the highest weight adjusted fraction (ensuring location constraint)

# INPUT: 
  # df: data frame containing the rejections made that should be filtered 
  # bigger_resolution: resolution number of the bigger of the two resolutions (integer)
  # smaller_resolution: resolution number of the smaller of the two resoulutions (integer)

# OUTPUT: 
  # final_df: returns for two given resolutions the dataframe containing the smallest group
return_selected_for_two_res <- function(df, bigger_resolution, smaller_resolution) {
  
  final_df <- c()
  
  groups_resolution_largest <- unique(df %>% filter(Resolution == bigger_resolution) %>% pull(Res_CHR_Group))
  
  if(length(groups_resolution_largest) > 0) {
    
    counter <- 0
    for(g in groups_resolution_largest) {
      
      
      # get SNPs in bigger resolution and fraction
      snps_in_largest_group <- df %>% filter(Resolution == bigger_resolution & Res_CHR_Group == g) %>% pull(SNP)
      frac_adjusted_largest_group <- unique(df %>% filter(Resolution == bigger_resolution & Res_CHR_Group == g) %>% pull(frac_adjusted))
      
      
      # now look at all groups and their fractions for the groups that the SNPs from the bigger resolution fall into in the smaller resolution
      # due to the <= here, we will keep those in lower level resolutions that we might have kept in previous iterations
      groups_in_smaller_resolution <- df %>% filter(Resolution <= smaller_resolution & SNP %in% snps_in_largest_group) %>% dplyr::select(Res_CHR_Group, frac_adjusted)
      groups_in_smaller_resolution <- unique(groups_in_smaller_resolution)
      
      # check if there is at least one element in the smaller groups that is bigger than the larger groups
      # if there is an element in the smaller group that (using the adjustment) is bigger than the larger group, keep the smaller group 
      # if the larger group has a larger adjusted fraction, keep the larger group 
      
      smaller_is_bigger <- (sum(groups_in_smaller_resolution$frac_adjusted >= frac_adjusted_largest_group) >= 1)
      
      # accounts for the case where there are no SNPs in lower-level resolution, adds the higher level one! (since sum() ... = 0 as groups_in_smaller_resolution$frac_adjusted = numeric(0) )
      # this could happen for example if the fraction in the lower lever was zero, then we will try to add the higher level one 
      # those cases where both the smaller level and the higher level are zero are not part of the dataset
      if(smaller_is_bigger) {
        #bigger_group <- groups_in_smaller_resolution[which.max(groups_in_smaller_resolution$frac_adjusted), ]$Res_CHR_Group
        #selected_res <- df %>% filter(Resolution <= smaller_resolution & Res_CHR_Group == bigger_group )
        selected_res <- df %>% filter(Resolution <= smaller_resolution & Res_CHR_Group %in% groups_in_smaller_resolution$Res_CHR_Group )
        final_df <- rbind(final_df, selected_res)
      }
      
      if(!smaller_is_bigger) {
        selected_res <- df %>% filter(Resolution == bigger_resolution & Res_CHR_Group == g)
        final_df <- rbind(final_df, selected_res)
      }
      
      
      
    }   
    
  }
  
  if(length(groups_resolution_largest) == 0) {
    final_df <- df
  }
  
  
  
  
  # returns for two given resolutions the dataframe containing the smallest group
  return(final_df)
  
}


# function to filter to outer nodes based on knockoff rejections 
# INPUT: 
  # df: data frame containing the rejections on individual level, including group identification 
  # resolutions: resolutions that have been used in the analysis 
  # group_identifiyer: variable name in df for the column with the group identification
# OUTPUT: 
  # outer nodes: data frame containing the outer nodes

filter_outer_node <- function(df, resolutions = all_resolutions, group_identifiyer = "Res_CHR_Group") {
  
  
  # for every resolution, look at all higher level rejections and keep only those with the smallest resolution
  
  outer_nodes <- c()
  
  for(r in resolutions) {
    
    
    #print(paste0("Outer Node Filtering, Resolution ", r))
    
    rej_groups_max_res <- df %>% filter(Resolution == r) %>% pull((!!as.name(group_identifiyer))) %>% unique()
    
    for(g in rej_groups_max_res) {
      
      # get SNPs in bigger resolution and fraction
      snps_in_largest_group <- df %>% filter((!!as.name(group_identifiyer)) == g) %>% pull(SNP)
      
      # find all groups containing these snps and keep all groups in the smallest resolution 
      groups_containing_snps_smallest_res <- df %>% filter(SNP %in% snps_in_largest_group) %>% 
        dplyr::select(Resolution,(!!as.name(group_identifiyer))) %>% unique() %>% filter(Resolution == min(Resolution))
      
      
      outer_nodes <- rbind(outer_nodes, groups_containing_snps_smallest_res)
      
    }
    
    
    
  }
  
  # there are duplicates, but only keep unique
  outer_nodes <- outer_nodes %>% unique()
  
  return(outer_nodes)
}





