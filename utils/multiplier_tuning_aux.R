
                  # MULTIPLIER TUNING #

# This file contains a function to generate a random matrix with fixed rowsum

# function to generate random matrix with fixed row-sum adding up to total_groups 
# INPUT: 
  # M: number of resolutions (integer)
  # total_groups: total number of groups (each row in output matrix sums up to total_groups, integer)
  # l: number of rows of matrix to generate 
  # add_extremes: indicator for whether to add extreme distributions (uniform distirbution, giving a single resolution all weight)
# OUTPUT: 
  # rounded_matrix: matrix of dimension (l + extremes) x M, where each row adds up to total_groups, filled with random integers
integer_matrix_fixed_rowsum_UKB <- function(M, total_groups, l = 100, add_extremes = TRUE, myseed = 2022) {
  
  # find matrix of allowed combinations 
  # this generates random (positive) real numbers, not integers, so round and then fix column sums
  set.seed(myseed)
  rounded_matrix <- t(round(Surrogate::RandVec(n = M, m = l, s = total_groups, a = 1, b = (total_groups) - (M-1))$RandVecOutput))
  
  rowsums <- rowSums(rounded_matrix)
  
  # make sure that all of them sum up to 1250
  for(i in 1:length(rowsums)) {
    
    if(rowsums[i] != total_groups) {
      
      # calculate difference 
      diff <- total_groups - rowsums[i]
      
      while(abs(diff) > 0) {
        # randomly choose a column and add or substract 1
        random_column <- sample.int(M, size = 1)
        if(diff < 0) {
          rounded_matrix[i, random_column] <- rounded_matrix[i, random_column] - 1
          rowsums[i] <- sum(rounded_matrix[i,] )
        } else {
          rounded_matrix[i, random_column] <- rounded_matrix[i, random_column] + 1
          rowsums[i] <- sum(rounded_matrix[i,] )
        }
        
        diff <- total_groups - rowsums[i]
        
      }
      
    }
    
  }
  
  # make sure that "extremes" are included 
  if(add_extremes) {
    all_equal_vector <- rep(rep(floor(total_groups/M), M))
    # randomly choose which one gets extra if not multiple
    randomly_sampled_group <- sample(M, 1)
    all_equal_vector[randomly_sampled_group] <- all_equal_vector[randomly_sampled_group] +  (total_groups - floor(total_groups/M)*M)
    
    floor(total_groups/M) + (total_groups - floor(total_groups/M)*M)
    diags <- diag(total_groups - (M-1), nrow = M, ncol = M)
    diags[upper.tri(diags)] <- 1
    diags[lower.tri(diags)] <- 1
    missing_to_total_groups <- (total_groups - floor(total_groups/M)*M)
    rounded_matrix <- rbind(rounded_matrix, 
                            all_equal_vector ,  # divide all equally,randomly choose which one gets extra if not multiple
                            diags)
    rounded_matrix <- unique(rounded_matrix)
  } 
  
  
  rounded_matrix <- unname(rounded_matrix)
  rounded_matrix <- unlist(rounded_matrix)
  return(rounded_matrix)
  
}
