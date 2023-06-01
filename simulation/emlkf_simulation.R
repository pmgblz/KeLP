#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)



#####################################################################
                           # EMLKF # 
# Runs simulation comparing e-MLKF to MLKF 
######################################################################

amp <- as.numeric(args[1])
if(is.na(amp)) amp <- 7

print(paste0("This is amp: ", amp))

alpha <- as.numeric(args[2])
if(is.na(alpha)) alpha <- 0.2


print(paste0("This is alpha: ", alpha))

# toggle whether to multiply the alpha used for the MLKF by the constant 
# to ensure that both the mlkf and the e-mlkf have the same FDR control threshold
alpha_ck_mlkf <- as.character(args[3])
alpha_ck_mlkf <- eval(parse(text = alpha_ck_mlkf)) 
if(is.na(alpha_ck_mlkf)) alpha_ck_mlkf <- FALSE


print(paste0("This is alpha_ck_mlkf: ", alpha_ck_mlkf))


ck <- 1.93
if(alpha_ck_mlkf) {
  alpha_mlkf <- alpha / ck
  alpha_emlkf <- alpha 
  ck_mult_ending = "mult"
} else {
  alpha_mlkf <- alpha 
  alpha_emlkf <- alpha
  ck_mult_ending = "not_mult"
}


print(paste0("This is alpha_mlkf: ", alpha_mlkf))
print(paste0("This is alpha_emlkf: ", alpha_emlkf))


gamma <- as.numeric(args[4])
if(is.na(gamma)) gamma <- alpha / 2

print(paste0("This is gamma: ", gamma))


nrep <- as.numeric(args[5])
if(is.na(nrep)) nrep <- 50

print(paste0("This is nrep: ", nrep))


nonzero_per_outcome <- as.numeric(args[6])  # how  many of the "largest" layer groups should be non-null 
if(is.na(nonzero_per_outcome)) nonzero_per_outcome <- 75
print(paste0("This is nonzero_per_outcome: ", nonzero_per_outcome))


number_non_null_big_groups <- as.numeric(args[7])  # now many non-zero elements should there be for each outcome 
if(is.na(number_non_null_big_groups)) number_non_null_big_groups <- 10
print(paste0("This is number_non_null_big_groups: ", number_non_null_big_groups))


block_sigma <- as.character(args[8])
block_sigma <- eval(parse(text = block_sigma)) 
if(is.na(block_sigma)) block_sigma <- FALSE
print(paste0("This is block_sigma: ", block_sigma))


tune_gamma_indicator <- as.character(args[9])
tune_gamma_indicator <- eval(parse(text = tune_gamma_indicator)) 
if(is.na(tune_gamma_indicator)) tune_gamma_indicator <- FALSE
print(paste0("This is tune_gamma: ", tune_gamma_indicator))


if(tune_gamma_indicator) {
  gamma_ending = "tuned"
} else {
  gamma_ending = "not_tuned"
}


tail_max <- as.character(args[10])
tail_max <- eval(parse(text = tail_max)) 
if(is.na(tail_max)) tail_max <- FALSE
print(paste0("This is tail_max: ", tail_max))


if(tail_max) {
  tail_max_ending = "max"
} else {
  tail_max_ending = "min"
}

same_gamma_each_m <- as.character(args[11])
same_gamma_each_m <- eval(parse(text = same_gamma_each_m)) 
if(is.na(same_gamma_each_m)) same_gamma_each_m <- FALSE
print(paste0("This is same_gamma_each_m: ", same_gamma_each_m))


if(same_gamma_each_m) {
  same_gamma_ending = "sm"
} else {
  same_gamma_ending = "dm"
}

######### 0. Baseline Problem Parameters ####### 

myseed <- 2022
Sys.sleep(1)  
set.seed(myseed)
Sys.sleep(1)  


n <- 4500 # number of samples 
p <- 2000 # number of predictors

Q = 1 # single outcome
overlap_pct = 1 # pct of non-zero snps that are identical across outcomes. since there is only a single outcome, it is always 1
FDP_hat_type = "kn+" # offset = 1

M = 2 # number of layers
alpha_emlkf <- rep(alpha_emlkf, M)
alpha_mlkf <- rep(alpha_mlkf, M)

group_size_with_individual <- c(1, 10) # group size

if(block_sigma) {
  rho_values <- c(0.1, 0.8) 
  pct_block_rho = 0.8
} else {
  rho_values = 0.3
}



knockoff_type = "fixed_equi"
# knockoff_method = "equi" # type of knockoff construction method FOR JULIA
statistic_type <- "group_LSM" # method to calcualte group importance statistics
offset = 1 # offset to use in knockoff 



############################## SET FUNCTIONS ####################


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


# LOAD ALL FUNCTIONS 
source(paste0(mydir, "utils/load_aux.R"))

# tell JuliaCall where is the Julia executable (the path can be found by typing `which julia` on terminal)
julia <- julia_setup(JULIA_HOME = my_julia_path)

# load Knockoffs.jl package
julia_library("Knockoffs")


########## 1. Generate Group Structure ########

set.seed(myseed)

group_struct_res <- generate_group_structure(group_size_with_individual, p, M) 

# groups is a data frame of dimension p x M 
# first column has group number in individual layer 
# second column has group number in second layer 
# ... and so on
groups <- group_struct_res[[1]]

# group_candidates_list is a list whith length == total number of groups (i.e. all individuals + all groups in other layers) 
# each list element contains the number of the individual elements that are part of each group
group_candidates_list <- group_struct_res[[2]]

group_data_frame <- group_struct_res[[3]]

######### 2. Generate Beta ############### 

# for multiple outcomes, based on overlap 
# use generate_beta_multiple_outcomes function 


# randomly samples number_non_null_big_groups
# for each sample non-null big group, randomly sample floor(overlap_pct*nonzero_per_outcome / number_non_null_big_groups) indices to be non-nulls 
# randomly choose one sampled non-null big_group to get the remainder 
# then randomly samples any index for any outcome that has not been sampled previously and belongs to the indices in the big groups
beta_res <- generate_beta_multiple_outcomes(p = p, Q = Q, M = M, groups = groups, 
                                        number_non_null_big_groups = number_non_null_big_groups,
                                        overlap_pct = overlap_pct,
                                        nonzero_per_outcome = nonzero_per_outcome, myseed = myseed) 

# this is a matrix with 0/1 indicators for whether beta is nonzero
beta_ind_nonzero <- beta_res$beta_ind_nonzero

# this is the final matrix also containing information on whether the effect is positive or negative
beta <- beta_res$beta_ind_random_sign


# non-null 
non_null_any <- (rowSums(beta_ind_nonzero) >= 1)
non_null_any_indices <- which(rowSums(beta_ind_nonzero) >= 1)

nonzero_groups_per_m = sapply(1:M, function(m)(unique(groups[non_null_any,m]))) 
print(nonzero_groups_per_m)

######### 3. Generate Sigma ############ 

# TUNE GAMMA 
# tune gamma if desired 
if(tune_gamma_indicator) {
  set.seed(myseed)
  final_gamma_combination <- emlkf_return_gamma(X = NA, Y = NA, 
                                                n = n, p = p, 
                                                amp = amp, 
                                                Sigma = Sigma, beta = beta,
                                                Q = Q, 
                                                groups = groups, knockoff_type = knockoff_type,
                                                alpha_emlkf = alpha_emlkf, 
                                                starting_gamma = alpha_emlkf[1] / 2, 
                                                FDP_hat_type = FDP_hat_type, 
                                                gamma_seq_params = c(gamma, gamma / 2, 0.05), 
                                                same_gamma_each_m = same_gamma_each_m, 
                                                tail_max = tail_max, 
                                                myseed = myseed, 
                                                generate_data = TRUE)
} else { 
  final_gamma_combination <- rep(gamma, M)
}

print("length / dimensions of gamma")
print(length(final_gamma_combination))


########## RUN SIMULATION ############

# create empty dataframe to store results
all_res <- c()
mycolnames <- c("Method", "Layer", "FDP", "Power", "Amplitude","Seed", "Gamma")

for(seedA in 1:nrep) {
  
  cat(sprintf("Running the %d-th rep.\n", seedA)) 
  
  
  # the actual data will be generated with a different seed
  seedB <- seedA + myseed
  set.seed(seedB)
  
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
                            myseed = seedB)
  
  Y <- as.matrix(data_res$Y)
  print(summary(Y))
  X <- as.matrix(data_res$X)
  print(summary(rowMeans(X)))
  Sigma <- as.matrix(data_res$Sigma)
  print(summary(rowMeans(Sigma)))

  # standardize columns of X
  X = scale(X, center = FALSE)/sqrt(n-1)

  # construct knockoffs at each layer without julia, asap group knockoffs, as in original paper 
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
  
  
  # RUN MULTILAYER KNOCKOFF FILTER
  cat("running mlkf ...")
  output_mlkf <- multilayer_knockoff_filter(X = X, Y = Y, W = W[[Q]], groups = groups, 
                                            alpha_mlkf = alpha_mlkf, 
                                                     FDP_hat_type = FDP_hat_type)
  cat("done.\n")
  
  # RUN E-MLKF
  cat("running e-mlkf ...")
  output_e_mlkf <- e_multilayer_knockoff_filter(X = X, Y = Y, W = W[[Q]], 
                                                groups = groups, 
                                                final_gamma_combination = final_gamma_combination, 
                                                alpha_emlkf = alpha_emlkf)
  cat("done.\n")
  

  # individual-layer rejection set MLKF
  S_hat_mlkf = output_mlkf$S_hat 
  S_hat_e_mlkf = output_e_mlkf$S_hat
  
  # rejection sets at each layer
  S_hats_mlkf = sapply(1:M, function(m)(unique(groups[S_hat_mlkf,m]))) 
  S_hats_e_mlkf = sapply(1:M, function(m)(unique(groups[S_hat_e_mlkf,m]))) 
  
  
  # get FDP and power
  FDP_mlkf = sapply(1:M, function(m)(length(setdiff(S_hats_mlkf[[m]], nonzero_groups_per_m[[m]]))/length(S_hats_mlkf[[m]])))
  power_mlkf = sapply(1:M, function(m)(length(intersect(S_hats_mlkf[[m]], nonzero_groups_per_m[[m]]))/length(nonzero_groups_per_m[[m]])))   
  
  FDP_e_mlkf = sapply(1:M, function(m)(length(setdiff(S_hats_e_mlkf[[m]], nonzero_groups_per_m[[m]]))/length(S_hats_e_mlkf[[m]])))
  power_e_mlkf = sapply(1:M, function(m)(length(intersect(S_hats_e_mlkf[[m]], nonzero_groups_per_m[[m]]))/length(nonzero_groups_per_m[[m]])))  
  
  # summarize in data frames
  layer = seq(1, M)
  seed = rep(seedA, M)
  
  method = rep("MLKF", M)
  summary_df_mlkf = data.frame(method, layer, FDP_mlkf, power_mlkf, amp, seed, gamma = rep(NA, M))
  colnames(summary_df_mlkf) <- mycolnames
  print("MLKF")
  print(summary_df_mlkf)
  
  method = rep("EMLKF", M)
  print(final_gamma_combination)
  print(dim(final_gamma_combination))
  summary_df_e_mlkf = data.frame(method, layer, FDP_e_mlkf, power_e_mlkf, amp, seed, gamma = final_gamma_combination)
  
 
  colnames(summary_df_e_mlkf) <- mycolnames
  print("EMLKF")
  print(summary_df_e_mlkf)
  
  all_res <- rbind(all_res, summary_df_mlkf, summary_df_e_mlkf)
  
}

print(all_res %>% group_by(Method, Layer) %>% mutate(fdp = ifelse(is.na(FDP), 0, FDP)) %>% 
        summarise(meanp = mean(Power), meanf = mean(fdp)))

## Save results
out_dir <- sprintf("%s/emlkf_res_tune_amp_%.1f_alpha_%.3f_nrep_%d_nnbg_%d_ck_%s_gamma_%s_ties_%s_gm_%s_2M_g_%s.csv", save_dir, 
                   amp, alpha,nrep, number_non_null_big_groups, ck_mult_ending, gamma_ending, tail_max_ending, same_gamma_ending, final_gammas)
write_csv(all_res, out_dir)
