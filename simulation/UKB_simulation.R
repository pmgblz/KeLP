#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)


######################################################################

# IMPORTANT: THE DATA IS NOT AVAILABLE. 
# The data is not ours. 
# We have applied for it from a source (the UK Biobank). 
# Interested parties can also apply to the same source.

######################################################################

# population can be one of whitenonbritish or british
# all samples are unrelated
my_population <- as.character(args[1])
if(is.na(my_population)) my_population <- "whitenonbritish"

amp <- as.numeric(args[2])
if(is.na(amp)) amp <- 7

alpha <- as.numeric(args[3])
if(is.na(alpha)) alpha <- 0.2

nrep <- as.numeric(args[4])
if(is.na(nrep)) nrep <- 20

Q <- as.numeric(args[5])
if(is.na(Q)) Q <- 1

sparsity <- as.numeric(args[6])  # how  many of the "largest" layer groups should be non-null 
if(is.na(sparsity)) sparsity <- 0.025

overlap_pct <- as.numeric(args[7]) # how many of the non-zero SNPs should be overlapping between multiple outcomes 
if(is.na(overlap_pct)) overlap_pct <- 1


partial_u <- as.character(args[8])
partial_u <- eval(parse(text = partial_u)) 
if(is.na(partial_u)) partial_u <- 1


tune_gamma_indicator <- as.numeric(args[9])
if(is.na(tune_gamma_indicator)) tune_gamma_indicator <- 0

cat(paste0("population is ", my_population, "\n"))
cat(paste0("amp is ", amp, "\n"))
cat(paste0("alpha is ", alpha, "\n"))
cat(paste0("nrep is ", nrep, "\n"))
cat(paste0("Q is ", Q, "\n"))
cat(paste0("sparsity is ", sparsity, "\n"))
cat(paste0("overlap_pct is ", overlap_pct, "\n"))
cat(paste0("partial u is ", partial_u, "\n"))

ncores = 1

print(paste0("Population: ", my_population, 
             "; N Cores: ", ncores))


# set your directory
if(dir.exists("")){
  mydir <- ""
  setwd(mydir)
  
  my_out_dir <- ""
}

## The directory to save the results
save_dir <- sprintf(paste0(mydir, "output/simulation"))
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}

# load all functions
source(paste0(mydir, "utils/load_aux.R"))


dfmax  <- 10000

# Set seed for cross-validated lasso statistics
myseed <- 2022
set.seed(myseed)

# partially using code from https://github.com/msesia/knockoffgwas/blob/master/knockoffgwas/utils/lasso.R

# Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(bigsnpr))
suppressMessages(library(knockoff))

chromosome_to_select <- 21


# --------------------------------------------------------------------------------------#
############  Load res 0 genotypes to define baseline group structure ########
# --------------------------------------------------------------------------------------#


cat("Finding baseline n, p")
baseline_resolution <- 0
bed_bim_fam_filename <- paste0("INSERT_FILEPATH_TO_DATA")
plinkfile = paste0(bed_bim_fam_filename,
                   ".bed")
if(!file.exists(paste0(bed_bim_fam_filename, ".bk"))) {
  x = snp_readBed2(plinkfile)
}

# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach(paste0(bed_bim_fam_filename,
                                ".rds"))

# Extract list of variants
map <- obj.bigSNP$map %>% as_tibble()
colnames(map) <- c("CHR", "SNP", "gd", "BP", "a1", "a2")
map <- map %>% dplyr::select(CHR, SNP, BP)


# Get aliases for useful slots
CHR_all <- obj.bigSNP$map$chromosome
CHR <- obj.bigSNP$map$chromosome[CHR_all == chromosome_to_select]
POS <- obj.bigSNP$map$physical.pos[CHR_all == chromosome_to_select]
SNP <- obj.bigSNP$map$marker.ID[CHR_all == chromosome_to_select]

# all genotypes
G_all <- obj.bigSNP$genotypes
# select only genotypes for specific chromosome
G <- big_copy(G_all, ind.col = which(CHR_all == chromosome_to_select))

# define n, p used for simulation
n <- nrow(G)
p <- (ncol(G) / 2) # divide by 2 because G contains the knockoffs as well

cat("done.\n")


# --------------------------------------------------------------------------------------#
############  Read genotypes for specific resolution ########
# --------------------------------------------------------------------------------------#


W.stats <- function(Z, knockoff) {
  importance <- abs(Z)
  z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
  zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
  w <- z-zk
}



simulate_y <- function(my_population, chromosome_to_select, beta_constructed, n, resolution = 0, myseed) {
  
  set.seed(myseed)
  
  cat(paste0("Processing genotypes resolution ", resolution, "..."))
  
  bed_bim_fam_filename <- paste0("INSERT_FILEPATH_TO_DATA")
  
  # load genotypes
  plinkfile = paste0(bed_bim_fam_filename,
                     ".bed")
  if(!file.exists(paste0(bed_bim_fam_filename, ".bk"))) {
    x = snp_readBed2(plinkfile)
  }
  
  
  # Attach the "bigSNP" object in R session
  obj.bigSNP <- snp_attach(paste0(bed_bim_fam_filename,
                                  ".rds"))
  
  # Extract list of variants
  map <- obj.bigSNP$map %>% as_tibble()
  colnames(map) <- c("CHR", "SNP", "gd", "BP", "a1", "a2")
  map <- map %>% dplyr::select(CHR, SNP, BP)
  
  # Extract list of subjects
  Subjects <- obj.bigSNP$fam %>% as_tibble()
  colnames(Subjects) <- c("FID", "IID", "X1", "X2", "sex", "X3")
  Subjects <- Subjects %>% dplyr::select(FID, IID) %>%
    mutate(FID=as.character(FID), IID=as.character(IID))
  
  # Get aliases for useful slots
  CHR_all <- obj.bigSNP$map$chromosome
  CHR <- obj.bigSNP$map$chromosome[CHR_all == chromosome_to_select]
  POS <- obj.bigSNP$map$physical.pos[CHR_all == chromosome_to_select]
  SNP <- obj.bigSNP$map$marker.ID[CHR_all == chromosome_to_select]
  
  # all genotypes
  G_all <- obj.bigSNP$genotypes
  # select only genotypes (includes knockoffs) for specific chromosome
  G <- big_copy(G_all, ind.col = which(CHR_all == chromosome_to_select))
  
  
  #####  ......Load variant partitions #######
  
  ## Load list of variants and partitions
  bim.file <- paste0(bed_bim_fam_filename, ".bim")
  Variants <- read.csv(bim.file, sep = "\t", header = FALSE)
  colnames(Variants) <-  c("CHR", "SNP", "X0", "BP", "X1", "X2")
  
  Variants <- Variants %>%
    mutate(Knockoff = ifelse(endsWith(SNP, ".k"), TRUE, FALSE), 
           SNP = sub(".k$", "", SNP)) %>% 
    dplyr::select(CHR, SNP, BP, Knockoff)
  
  # Import Group Information 
  myfilename <- paste0("INSERT_FILEPATH_TO_DATA")
  grp.file.chr <- read.csv(myfilename, sep="") 
  grp.file.chr$CHR <- chromosome_to_select
  
  Variants <- Variants %>% 
    mutate(CHR = as.integer(CHR)) %>%
    dplyr::filter(CHR == chromosome_to_select) %>%
    left_join(grp.file.chr, by = c("SNP", "CHR"))
  
  # get indices of original X (not knockoff)
  indices_actual_genotypes <- which(Variants$Knockoff == FALSE)
  
  group_bp_min_max <- Variants %>% 
    group_by(CHR, Group) %>% 
    summarise(BP.min = min(BP), BP.max = max(BP))
  
  # Compute scaling factor for the genotypes
  scaler <- big_scale()
  G.scale <- scaler(G)
  scaling.factors <- G.scale$scale
  
  # select only true genotypes (not the knockoffs) for specific chromosome
  # variants is already filtered to the chromosome of interest 
  # so also use the filtered genotypes here
  X <- big_copy(G, ind.col = indices_actual_genotypes)
  
  # get covariates
  # Load phenotype table
  covariates_file <- read.csv(paste0("INSERT_FILEPATH_TO_DATA"))
  
  # extract covariates 
  Covariates <- covariates_file %>% dplyr::select(age, age_squared, sex, PC1, PC2, PC3, PC4, PC5)
  Covariates <- as.matrix(Covariates)
  
  
  ###### ..... SIMULATE PHENOTYPE ########## 
  
  
  Y <- matrix(NA, nrow = n, ncol = Q)
  
  # generate data if using signal amplitude
  # beta was fixed before, amp / sqrt(n) is nonrandom
  for(q in 1:Q) {
    beta_true  <- as.matrix(amp * beta_constructed[, q] / sqrt(n), nrow = n, ncol = 1)
    Y[, q] <- big_prodMat(X, beta_true) + rnorm(n)
  }
  
  my_snr <- var(big_prodMat(X, beta_true))
  
  return(list(X = X, Y = Y, my_snr = my_snr))
  
}



# function to load genotypes and run lasso for particular resolution
lasso_imp_stats_specific_resolution <- function(Y, my_population, chromosome_to_select, beta_constructed, gamma, alpha, Q, n, p, resolution, myseed) {
  
  set.seed(myseed)
  
  cat(paste0("Processing genotypes resolution ", resolution, "..."))
  
  bed_bim_fam_filename <- paste0("INSERT_FILEPATH_TO_DATA")
  
  # load genotypes
  plinkfile = paste0(bed_bim_fam_filename,
                     ".bed")
  if(!file.exists(paste0(bed_bim_fam_filename, ".bk"))) {
    x = snp_readBed2(plinkfile)
  }
  
  
  # Attach the "bigSNP" object in R session
  obj.bigSNP <- snp_attach(paste0(bed_bim_fam_filename,
                                  ".rds"))
  
  # Extract list of variants
  map <- obj.bigSNP$map %>% as_tibble()
  colnames(map) <- c("CHR", "SNP", "gd", "BP", "a1", "a2")
  map <- map %>% dplyr::select(CHR, SNP, BP)
  
  # Extract list of subjects
  Subjects <- obj.bigSNP$fam %>% as_tibble()
  colnames(Subjects) <- c("FID", "IID", "X1", "X2", "sex", "X3")
  Subjects <- Subjects %>% dplyr::select(FID, IID) %>%
    mutate(FID=as.character(FID), IID=as.character(IID))
  
  # Get aliases for useful slots
  CHR_all <- obj.bigSNP$map$chromosome
  CHR <- obj.bigSNP$map$chromosome[CHR_all == chromosome_to_select]
  POS <- obj.bigSNP$map$physical.pos[CHR_all == chromosome_to_select]
  SNP <- obj.bigSNP$map$marker.ID[CHR_all == chromosome_to_select]
  
  # all genotypes
  G_all <- obj.bigSNP$genotypes
  # select only genotypes (includes knockoffs) for specific chromosome
  G <- big_copy(G_all, ind.col = which(CHR_all == chromosome_to_select))
  
  
  #####  ......Load variant partitions #######
  
  ## Load list of variants and partitions
  bim.file <- paste0(bed_bim_fam_filename, ".bim")
  Variants <- read.csv(bim.file, sep = "\t", header = FALSE)
  colnames(Variants) <-  c("CHR", "SNP", "X0", "BP", "X1", "X2")
  
  Variants <- Variants %>%
    mutate(Knockoff = ifelse(endsWith(SNP, ".k"), TRUE, FALSE), 
           SNP = sub(".k$", "", SNP)) %>% 
    dplyr::select(CHR, SNP, BP, Knockoff)
  
  # Import Group Information 
  myfilename <- paste0("INSERT_FILEPATH_TO_DATA")
  grp.file.chr <- read.csv(myfilename, sep="") 
  grp.file.chr$CHR <- chromosome_to_select
  
  Variants <- Variants %>% 
    mutate(CHR = as.integer(CHR)) %>%
    dplyr::filter(CHR == chromosome_to_select) %>%
    left_join(grp.file.chr, by = c("SNP", "CHR"))
  
  # get indices of original X (not knockoff)
  indices_actual_genotypes <- which(Variants$Knockoff == FALSE)
  
  group_bp_min_max <- Variants %>% 
    group_by(CHR, Group) %>% 
    summarise(BP.min = min(BP), BP.max = max(BP))
  
  # Compute scaling factor for the genotypes
  scaler <- big_scale()
  G.scale <- scaler(G)
  scaling.factors <- G.scale$scale
  
  # get covariates
  # Load phenotype table
  covariates_file <- read.csv(paste0("INSERT_FILEPATH_TO_DATA"))
  
  # extract covariates 
  Covariates <- covariates_file %>% dplyr::select(age, age_squared, sex, PC1, PC2, PC3, PC4, PC5)
  Covariates <- as.matrix(Covariates)
  
  # fit lasso for each outcome
  W <- list()
  rejected_knockoff <- list()
  taus_evals <- list()
  taus_knockoff <- list()
  fracs <- matrix(NA, nrow = length(unique(Variants$Group)) , ncol = Q)
  fracs_knockoff <- matrix(NA, nrow = length(unique(Variants$Group)) , ncol = Q)
  
  
  cat(paste0("Running the lasso for resolution ", resolution))
  
  set.seed(myseed)
  
  for(q in 1:Q) {
    
    lasso.fit <- big_spLinReg(G, 
                              y.train=Y[, q],
                              covar.train=Covariates,
                              pf.covar = rep(0, ncol(Covariates)), # don't penalize covariates
                              dfmax=dfmax, 
                              ncores=ncores)
    
    
    # Extract beta from each fold and combine them
    beta <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)
    
    # Separate the coefficients of the genetic variants from the coefficients of the covariates
    beta.variants <- beta[1:ncol(G),]
    beta.covariates <- beta[(ncol(G)+1):nrow(beta),]
    
    # Undo scaling of lasso coefficients
    beta.variants <- beta.variants * scaling.factors
    Beta <- cbind(tibble(CHR=Variants$CHR,
                         SNP=Variants$SNP, BP=Variants$BP, Knockoff=Variants$Knockoff),
                  as_tibble(beta.variants)) %>% as_tibble()
    colnames(Beta) <- c("CHR", "SNP", "BP", "Knockoff", paste("K", seq(ncol(beta.variants)),sep=""))
    Beta <- Beta %>%
      mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
             Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
      dplyr::select(CHR, SNP, BP, Knockoff, Z)
    
    # Extract the estimated coefficients
    Lasso.res <- Beta %>%
      inner_join(Variants, by = c("CHR", "SNP", "BP", "Knockoff")) %>%
      dplyr::select(CHR, SNP, BP, Z, Group, Knockoff)
    
    
    # Compute the knockoff statistics
    Stats <- Lasso.res %>%
      dplyr::select("CHR", "Group", "SNP", "BP", "Knockoff", "Z") %>%
      group_by(CHR, Group) %>%
      summarize(W = W.stats(abs(Z),Knockoff),
                Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
                Size=n()) %>%
      ungroup() %>%
      dplyr::select(CHR, Group, SNP.lead, BP.lead, Size, W) %>% 
      arrange((abs(Group))) 
    
    # add in BP min and BP max information
    Stats <- left_join(Stats, group_bp_min_max, by = c("CHR", "Group"))
    
    # save importance statistics
    W[[q]] <- Stats$W
    
    taus_evals[[q]] = knockoff.threshold(W[[q]], fdr=gamma, offset=1)
    
    # calculate e-values
    fracs[, q] <- ((W[[q]] >=  taus_evals[[q]] ) / (1 + sum(W[[q]] <= -taus_evals[[q]] )))
    
    # rejected knockoff
    taus_knockoff[[q]] = knockoff.threshold(W[[q]], fdr=alpha, offset=1)
    rejected_knockoff[[q]] <- ifelse(W[[q]] >= taus_knockoff[[q]], 1, 0)
    fracs_knockoff[, q] <- ((W[[q]] >=  taus_knockoff[[q]] ) / (1 + sum(W[[q]] <= -taus_knockoff[[q]] )))
  }
  
  
  cat("done\n")
  
  return(list(W = W, fracs = fracs, fracs_knockoff = fracs_knockoff, rejected_knockoff = rejected_knockoff))
  
  
}




# --------------------------------------------------------------------------------------#
############  select non-zero groups ########
# --------------------------------------------------------------------------------------#


# create group file for all resolutions 
# Import Group Information for resolution 0

myfilename <- paste0("INSERT_FILEPATH_TO_DATA")
grp.file <- read.csv(myfilename, sep="") 
colnames(grp.file) <- c("SNP", paste0("Group_res", 0))

remaining_resolutions <- seq(1, 6)
for(r in remaining_resolutions) {
  
  myfilename <- paste0("INSERT_FILEPATH_TO_DATA")
  myfile <- read.csv(myfilename, sep="") 
  colnames(myfile) <- c("SNP", paste0("Group_res", r))
  
  grp.file <- left_join(grp.file, myfile, by = "SNP")
  
}


groups <- grp.file[, 2:8]
considered_resolutions <- c(0, 2, 4, 6)
groups <- groups[, 1 + considered_resolutions]


# CREATE BETA
beta_res <- random_beta_generation(p = p, Q = Q, sparsity = sparsity, overlap_pct = overlap_pct, myseed = myseed)

beta_constructed <- beta_res$beta_ind_random_sign
nonzero_matrix <- beta_res$nonzero_matrix
nonzero_per_q <- beta_res$nonzero_per_q
nonzero_any_indiv <- beta_res$nonzero_any_indiv
print("NUMBER NONZERO")
print(length(nonzero_any_indiv))

# nonzero per layer
nonzero_groups_per_layer <-sapply(1:length(considered_resolutions), function(m)(unique(groups[nonzero_any_indiv,m]))) 

# nonzero snps
nonzero_snps <- grp.file[grp.file$Group_res0 %in% nonzero_groups_per_layer[[1]], ]$SNP

# get vector containing nonzero res group identifiers
nonzero_res_group_per_layer <- c() 
counter = 0
nonzero_res_group_per_layer_list <- list()
for(m in considered_resolutions) {
  counter = counter + 1
  nonzero_res_group_per_layer_list[[counter]] <- paste0("res_", m, "_group_", nonzero_groups_per_layer[[counter]])
  nonzero_res_group_per_layer  <- c(nonzero_res_group_per_layer, paste0("res_", m, "_group_", nonzero_groups_per_layer[[counter]]))
}

# get snp-level file for considered resolutions only
grp.file.considered.res <- grp.file %>% 
  dplyr::select("SNP", ends_with(paste0("_res", considered_resolutions)))


cand_groups <- list() 
group_data_frame <- c()
inverse_group_size_by_res <- list()

for(i in 1:length(considered_resolutions)) {
  
  groups_specific_resolution <- grp.file.considered.res %>% 
    dplyr::select("SNP", ends_with(paste0("_res", considered_resolutions[i]))) %>% 
    mutate(Resolution = considered_resolutions[i])
  
  colnames(groups_specific_resolution) <- c("SNP", "group_number", "Resolution")
  
  groups_specific_resolution <- groups_specific_resolution %>% 
    mutate(group = paste0("group_", group_number)) %>% 
    mutate(Res_Group = paste0("res_", Resolution, "_", group)) %>% 
    mutate(res_Group = Res_Group)
  
  group_data_frame <- rbind(group_data_frame, groups_specific_resolution)
  
  
  list_specific_resolution <- split(groups_specific_resolution, groups_specific_resolution$group)
  list_specific_resolution <- lapply(list_specific_resolution, function(x) create_list_element(list_element = x))
  list_specific_resolution <- unname(list_specific_resolution)
  
  # order by group number
  list_specific_resolution <-  list_specific_resolution[ order ( sapply(list_specific_resolution, "[[", "group_number") )]
  
  # save group size
  inverse_group_size_by_res[[i]] <- sapply(list_specific_resolution, inverse_size)
  
  cand_groups <- c(cand_groups, list_specific_resolution)
  
}

# get group data frame info (snp level information with group and res_group)
group_data_frame_nonzero <- group_data_frame %>% dplyr::filter(Res_Group %in% nonzero_res_group_per_layer)
total_number_groups <- length(cand_groups)

# set gammas to be alpha / 4
gammas <- c(0.05, 0.05, 0.05, 0.05)


# --------------------------------------------------------------------------------------#
############  RUN SIMULATION ########
# --------------------------------------------------------------------------------------#


all_res_by_resolutions <- c()
all_res_across_resolutions <- c()
all_res_by_q <- c()
my_snr <- c()

## Run both procedures for multiple runs conditional on the realized (X,Y)
for(seedA in 1:nrep){
  
  cat(sprintf("Running the %d-th rep.\n", seedA)) 
  
  seedB <- seedA + myseed
  set.seed(seedB)
  
  
  # SIMULATE Y
  sim_y_res <- simulate_y(my_population = my_population, chromosome_to_select = chromosome_to_select, 
                          beta_constructed = beta_constructed,n, resolution = 0, myseed = seedB) 
  
  set.seed(seedB)
  
  Y <- sim_y_res$Y
  
  
  
  W <- list()
  fracs <- list()
  rejected_knockoff <- list()
  fracs_knockoff <- list()
  
  # run optimization problem
  old_time <- Sys.time() # get start time
  
  alpha = 0.2
  
  for(i in 1:length(considered_resolutions)) {
    
    
    set.seed(seedB)
    # this simulates Y, runs the lasso, and gives W as output
    res_i <- lasso_imp_stats_specific_resolution(Y = Y, 
                                                 my_population = my_population, chromosome_to_select = chromosome_to_select, 
                                                 beta_constructed = beta_constructed, gamma = gammas[i], alpha = 0.2, Q = 1, n = n, p = p,
                                                 resolution = considered_resolutions[i], myseed = seedB)
    
    W[[i]] <- res_i$W
    fracs[[i]] <- res_i$fracs
    rejected_knockoff[[i]] <- res_i$rejected_knockoff
    fracs_knockoff[[i]] <- res_i$fracs_knockoff
    
    
  }
  
  new_time <- Sys.time() - old_time # calculate difference
  print("RUNTIME LASSOs:")
  print(new_time) 
  
  
  ############### ... Power for Knockoff Separately at each Resolution ##################
  # power for knockoff at each resolution separately
  knockoff_by_resolution_rejections <- lapply(rejected_knockoff, function(x) which(x[[1]] == 1) - 1) # groups start at 0, so substract 1
  
  if(length(knockoff_by_resolution_rejections) > 0) {
    
    
    knockoff_by_resolution_fdr =  matrix(NA, nrow = 1, ncol = length(considered_resolutions)) 
    knockoff_by_resolution_true_power =  matrix(NA, nrow = 1, ncol = length(considered_resolutions)) 
    knockoff_by_resolution_power =  matrix(NA, nrow = 1, ncol = length(considered_resolutions)) 
    knockoff_by_resolution_snps_implicated =  matrix(NA, nrow = 1, ncol = length(considered_resolutions)) 
    knockoff_by_resolution_frac_correct =  matrix(NA, nrow = 1, ncol = length(considered_resolutions)) 
    
    counter = 0
    
    for(m in considered_resolutions) {
      
      counter = counter + 1
      
      if(length(knockoff_by_resolution_rejections[[counter]]) > 0) {
        
        nonzero_specific_layer <-  nonzero_res_group_per_layer_list[[counter]]
        knockoff_rejected_res_group_per_layer  <- paste0("res_", m, "_group_", knockoff_by_resolution_rejections[[counter]])
        
        knockoff_by_resolution_fdr[, counter] = length(setdiff(knockoff_rejected_res_group_per_layer, nonzero_specific_layer)) / length(knockoff_rejected_res_group_per_layer)
        knockoff_by_resolution_true_power[, counter] =  length(intersect(knockoff_rejected_res_group_per_layer, nonzero_specific_layer)) / length(nonzero_specific_layer)
        
        
        knockoff_true_group_rejections <- intersect(knockoff_rejected_res_group_per_layer, nonzero_specific_layer)
        # get the size for each of these true rejections
        group_data_frame_nonzero %>% 
          dplyr::filter(Res_Group %in% knockoff_true_group_rejections) %>% 
          group_by(Res_Group) %>% 
          count() %>% 
          mutate(inverse_size = 1 / n) -> inverse_size_true_knockoff_rejections
        
        # ALL REJECTIONS
        group_data_frame %>% 
          dplyr::filter(Res_Group %in% knockoff_rejected_res_group_per_layer) %>% 
          pull(SNP) %>% 
          unique() -> knockoff_rejected_snps
        
        
        # sum of inverse size of ture rejections divided by total number of nonzero snps (only counts the true rejections!)
        knockoff_by_resolution_power[, counter] <- sum(inverse_size_true_knockoff_rejections$inverse_size) / length(nonzero_groups_per_layer[[1]])
        
        knockoff_by_resolution_snps_implicated[, counter] <- length(knockoff_rejected_snps)
        
        knockoff_by_resolution_frac_correct[, counter] = length(intersect(knockoff_rejected_snps, nonzero_snps)) / length(nonzero_groups_per_layer[[1]])
        
        
      }
      
    }
    
    
  } else {
    knockoff_by_resolution_fdr <- rep(0, length(considered_resolutions))
    knockoff_by_resolution_power <- rep(0, length(considered_resolutions))
    knockoff_by_resolution_true_power <- rep(0, length(considered_resolutions))
    knockoff_by_resolution_snps_implicated <- rep(0, length(considered_resolutions))
    knockoff_by_resolution_frac_correct <- rep(0, length(considered_resolutions))
  }
  
  
  
  ############## ....Knockoff outer nodes #############
  
  # reformat knockoff rejections to pass into outer node filter
  rej_knockoff_res <- c()
  
  counter = 0
  for(m in considered_resolutions) {
    
    
    counter = counter + 1
    rejected_knockoff_group_df_res <- data.frame(Group = sort(unique(grp.file.considered.res[, 1 + counter])), 
                                                 rejected_knockoff = rejected_knockoff[[counter]][[1]], 
                                                 Resolution = considered_resolutions[counter]) %>%
      mutate(Res_Group = paste0("res_", Resolution, "_group_", Group))
    
    
    snp_group <- grp.file.considered.res[, c(1, 1 + counter)]
    colnames(snp_group) <- c("SNP", "Group")
    
    snp_group <- left_join(snp_group, rejected_knockoff_group_df_res, by = "Group")
    
    rej_knockoff_res <- rbind(rej_knockoff_res, snp_group)
    
    
  }
  
  
  rej_knockoff_res <- rej_knockoff_res %>% dplyr::filter(rejected_knockoff == 1)
  
  # filter to outer nodes
  outer_nodes_knockoff <- filter_outer_node(rej_knockoff_res, resolutions = considered_resolutions, group_identifiyer = "Res_Group")
  
  
  ################## ............... Power knockoff outer ######################
  # calculate power, fdr, snps implicated, frac correct
  selected_outer_knockoff_groups <- outer_nodes_knockoff$Res_Group
  
  if(length(selected_outer_knockoff_groups) > 0) {
    # get power, fdr and frac correct rejections
    knockoff_outer_fdr = length(setdiff(selected_outer_knockoff_groups, nonzero_res_group_per_layer)) / length(selected_outer_knockoff_groups)
    
    
    true_group_rejections <- intersect(selected_outer_knockoff_groups, nonzero_res_group_per_layer)
    # get the size for each of these true rejections
    group_data_frame_nonzero %>% 
      dplyr::filter(Res_Group %in% true_group_rejections) %>% 
      group_by(Res_Group) %>% 
      count() %>% 
      mutate(inverse_size = 1 / n) -> inverse_size_true_knockoff_outer_rejections
    
    group_data_frame %>% 
      dplyr::filter(Res_Group %in% selected_outer_knockoff_groups) %>% 
      pull(SNP) %>% 
      unique() -> knockoff_outer_rejected_snps
    
    group_data_frame %>% 
      dplyr::filter(Res_Group %in% selected_outer_knockoff_groups) %>% 
      group_by(Res_Group) %>%
      dplyr::summarise(group_size = n()) -> knockoff_outer_group_size_rejected
    
    
    knockoff_outer_avg_size_rejected_group = mean(knockoff_outer_group_size_rejected$group_size)
    knockoff_outer_median_size_rejected_group = median(knockoff_outer_group_size_rejected$group_size)
    
    
    # sum of inverse size of ture rejections divided by total number of nonzero snps
    knockoff_outer_power <- sum(inverse_size_true_knockoff_outer_rejections$inverse_size) / length(nonzero_groups_per_layer[[1]])
    
    knockoff_outer_snps_implicated <- length(knockoff_outer_rejected_snps)
    
    knockoff_outer_frac_correct = length(intersect(knockoff_outer_rejected_snps, nonzero_snps)) / length(nonzero_groups_per_layer[[1]])
    
    ############# ................e-BH on knockoff / global outer nodes #########
    
    # get fractions for each 
    
    fractions_knockoff_outer_rejected_per_layer <- c()
    
    for(i in 1:length(considered_resolutions)) {
      
      # unique groups, sorted
      group_data_frame %>% 
        dplyr::filter(Resolution == considered_resolutions[i]) %>%
        pull(group_number) %>% 
        unique() %>% 
        sort() -> groups_res_sorted
      
      rejected_frac_df <- data.frame(group = groups_res_sorted, 
                                     Res_Group = paste0("res_", considered_resolutions[i], "_group_", groups_res_sorted),
                                     fracs = fracs_knockoff[[i]], 
                                     eval = fracs_knockoff[[i]] * length(groups_res_sorted)) %>% 
        dplyr::filter(Res_Group %in% selected_outer_knockoff_groups) %>% 
        dplyr::select(Res_Group, fracs, eval)
      
      fractions_knockoff_outer_rejected_per_layer <- rbind(fractions_knockoff_outer_rejected_per_layer, rejected_frac_df)
      
      
    }
    
    # now run e-BH on these
    
    # Q = 1, so it will always be u = 1 and alpha not adjusted for partial conjunction (only for selection)
    number_groups_outer <- length(unique(fractions_knockoff_outer_rejected_per_layer$Res_Group))
    alpha_mod <- alpha * number_groups_outer / total_number_groups
    
    ebh_outer_nodes_rej_indices <- run_eBH(fractions_knockoff_outer_rejected_per_layer$eval, alpha_mod, 
                                           number_groups = number_groups_outer, use_fractions = FALSE)
    
    ######################## ............... Power for eBH on knockoff outer #################
    
    if(!is.na(ebh_outer_nodes_rej_indices[1])) {
      
      ebh_knockoff_outer_rejected_res_group <- fractions_knockoff_outer_rejected_per_layer[ebh_outer_nodes_rej_indices,]$Res_Group
      
      # get power, fdr and frac correct rejections
      ebh_knockoff_outer_fdr = length(setdiff(ebh_knockoff_outer_rejected_res_group, nonzero_res_group_per_layer)) / length(ebh_knockoff_outer_rejected_res_group)
      
      
      true_group_rejections_knockoff_ebh <- intersect(ebh_knockoff_outer_rejected_res_group, nonzero_res_group_per_layer)
      # get the size for each of these true rejections
      group_data_frame_nonzero %>% 
        dplyr::filter(Res_Group %in% true_group_rejections_knockoff_ebh) %>% 
        group_by(Res_Group) %>% 
        count() %>% 
        mutate(inverse_size = 1 / n) -> inverse_size_true_ebh_knockoff_outer_rejections
      
      group_data_frame %>% 
        dplyr::filter(Res_Group %in% ebh_knockoff_outer_rejected_res_group) %>% 
        pull(SNP) %>% 
        unique() -> ebh_knockoff_outer_rejected_snps
      
      # sum of inverse size of ture rejections divided by total number of nonzero snps
      ebh_knockoff_outer_power <- sum(inverse_size_true_ebh_knockoff_outer_rejections$inverse_size) / length(nonzero_groups_per_layer[[1]])
      
      ebh_knockoff_outer_snps_implicated <- length(ebh_knockoff_outer_rejected_snps)
      
      ebh_knockoff_outer_frac_correct = length(intersect(ebh_knockoff_outer_rejected_snps, nonzero_snps)) / length(nonzero_groups_per_layer[[1]])
      
      
    } else {
      ebh_knockoff_outer_fdr = 0
      ebh_knockoff_outer_power = 0
      ebh_knockoff_outer_snps_implicated = 0
      ebh_knockoff_outer_frac_correct = 0
    }
    
    
  } else {
    knockoff_outer_fdr = 0
    knockoff_outer_power = 0
    knockoff_outer_frac_correct = 0
    knockoff_outer_snps_implicated = 0
    knockoff_outer_avg_size_rejected_group = 0
    knockoff_outer_median_size_rejected_group = 0
    ebh_knockoff_outer_fdr = 0
    ebh_knockoff_outer_power = 0
    ebh_knockoff_outer_snps_implicated = 0
    ebh_knockoff_outer_frac_correct = 0
  }
  
  ############# ....Focused e-BH / KELP ###############
  
  # run focused e-BH (same as kelp in this case)
  cand_groups_kelp <- cand_groups
  
  for(g in 1:length(cand_groups_kelp)) {
    cand_groups_kelp[[g]]$fracs <- unlist(fracs)[g]
  }
  
  # filter out groups with zero e-values (or fractions)
  cand_groups_kelp <- cand_groups_kelp[sapply(cand_groups_kelp, function(x) x$fracs > 0)]
  
  # get individual "snps" in the non-zero groups
  individual_snps_candidates <-  unique(unlist(sapply(cand_groups_kelp, `[[`, "SNP")))
  
  # number groups in optimization procedure
  n_groups_filtered <- length(cand_groups_kelp)
  
  # number of unique snps 
  n_individual_hypotheses_filtered <- length(individual_snps_candidates)
  
  # run optimization problem
  
  if(n_groups_filtered > 0) {
    cat("Running kelp ...")
    res <- kelp(group_candidates_list = cand_groups_kelp, 
                individual_snps_candidates = individual_snps_candidates, 
                single_u = TRUE,
                alpha = alpha, 
                verbose = F,
                weighted = T, 
                fractions_only = TRUE,
                M = length(considered_resolutions), 
                n_individual_hypotheses_filtered = n_individual_hypotheses_filtered, 
                n_groups_filtered = n_groups_filtered, 
                total_groups = length(cand_groups), 
                Q = Q,
                u_seq = 1) 
    cat("done\n.")
    kelp_detections <- res[[1]]
    
    # get selected groups
    selected_kelp_groups <- sapply(kelp_detections, `[[`, "Res_Group")
  } else {
    selected_kelp_groups <- NULL
  }
  
  
  if(length(selected_kelp_groups) > 0) {
    
    kelp_rejected_res_group <- selected_kelp_groups
    
    # get power, fdr and frac correct rejections
    kelp_fdr = length(setdiff(kelp_rejected_res_group, nonzero_res_group_per_layer)) / length(kelp_rejected_res_group)
    
    
    true_group_rejections_kelp <- intersect(kelp_rejected_res_group, nonzero_res_group_per_layer)
    # get the size for each of these true rejections
    group_data_frame_nonzero %>% 
      dplyr::filter(Res_Group %in% true_group_rejections_kelp) %>% 
      group_by(Res_Group) %>% 
      count() %>% 
      mutate(inverse_size = 1 / n) -> inverse_size_true_kelp_rejections
    
    group_data_frame %>% 
      dplyr::filter(Res_Group %in% kelp_rejected_res_group) %>% 
      pull(SNP) %>% 
      unique() -> kelp_rejected_snps
    
    group_data_frame %>% 
      dplyr::filter(Res_Group %in% kelp_rejected_res_group) %>% 
      dplyr::group_by(Res_Group) %>% 
      dplyr::summarise(group_size = n()) -> group_size_kelp
    
    kelp_outer_avg_size_rejected_group <- mean(group_size_kelp$group_size)
    kelp_median_size_rejected_group <- median(group_size_kelp$group_size)
    
    
    # sum of inverse size of ture rejections divided by total number of nonzero snps
    kelp_power <- sum(inverse_size_true_kelp_rejections$inverse_size) / length(nonzero_groups_per_layer[[1]])
    
    kelp_snps_implicated <- length(kelp_rejected_snps)
    
    kelp_frac_correct = length(intersect(kelp_rejected_snps, nonzero_snps)) / length(nonzero_groups_per_layer[[1]])
    
    
    
    
  } else { 
    kelp_fdr <- 0
    kelp_power <- 0
    kelp_snps_implicated <- 0
    kelp_frac_correct <- 0
    kelp_outer_avg_size_rejected_group = 0
    kelp_median_size_rejected_group = 0
  }
  
  
  
  
  
  
  
  # by Resolution 
  summary_by_resolution = data.frame(level = rep("by_resolution", length(considered_resolutions)),
                                     method = c(rep("knockoff", length(considered_resolutions))),
                                     Resolution = considered_resolutions,
                                     fdp = c(knockoff_by_resolution_fdr),
                                     power = c(knockoff_by_resolution_power),
                                     power_true = c(knockoff_by_resolution_true_power),
                                     snps_implicated = c(knockoff_by_resolution_snps_implicated), 
                                     frac_correct_rejections = c(knockoff_by_resolution_frac_correct),
                                     amp = rep(amp, length(considered_resolutions)), 
                                     seed = rep(seedA, length(considered_resolutions)), 
                                     Q = rep(Q, length(considered_resolutions)), 
                                     sparsity = rep(sparsity, length(considered_resolutions)), 
                                     overlap_pct = rep(overlap_pct, length(considered_resolutions)), 
                                     alpha = rep(alpha, length(considered_resolutions)), 
                                     gamma = gammas)
  
  
  summary_across_resolutions = data.frame(level = "across_resolution",
                                          method = c("kelp", "knockoffouter", "ebhknockout"),
                                          fdp = c(kelp_fdr,knockoff_outer_fdr, ebh_knockoff_outer_fdr),
                                          power = c(kelp_power, knockoff_outer_power, ebh_knockoff_outer_power),
                                          snps_implicated = c(kelp_snps_implicated, knockoff_outer_snps_implicated, ebh_knockoff_outer_snps_implicated),
                                          frac_correct_rejections = c(kelp_frac_correct,knockoff_outer_frac_correct, ebh_knockoff_outer_frac_correct),
                                          avg_group_size_rejected = c(kelp_outer_avg_size_rejected_group, knockoff_outer_avg_size_rejected_group, NA), 
                                          median_group_size_rejected = c(kelp_median_size_rejected_group, knockoff_outer_median_size_rejected_group, NA),
                                          amp = amp, 
                                          seed = seedA, 
                                          Q = Q, 
                                          sparsity = sparsity, 
                                          overlap_pct = overlap_pct, 
                                          alpha = alpha, 
                                          gamma = paste(gammas, collapse = "_"))
  
  
  snr_df <- data.frame(my_snr = sim_y_res$my_snr, 
                       amp = amp, 
                       seed = seedA, 
                       Q = Q, 
                       sparsity = sparsity, 
                       overlap_pct = overlap_pct, 
                       alpha = alpha)
  
  my_snr <- rbind(my_snr, snr_df)
  
  
  all_res_by_resolutions <- rbind(all_res_by_resolutions, summary_by_resolution)
  all_res_across_resolutions <- rbind(all_res_across_resolutions, summary_across_resolutions)
  print(all_res_across_resolutions)
  
}  

## Save the outcomes
out_dir <- sprintf("%s/UKB_kelp_res_by_M_%d_Q_%d_amp_%.1f_nrep_%d_ol_%.1f_g_%s_sp%.4f_al%.4f_u%d_c%d_tg%d.csv", save_dir, length(considered_resolutions), Q,
                   amp, nrep, overlap_pct, paste(gammas, collapse = "_"), sparsity, alpha, partial_u, chromosome_to_select, tune_gamma_indicator)
write_csv(data.frame(all_res_by_resolutions), out_dir)

out_dir <- sprintf("%s/UKB_kelp_res_across_M_%d_Q_%d_amp_%.1f_nrep_%d_ol_%.1f_g_%s_sp%.4f_al%.4f_u%d_c%d_tg%d.csv", save_dir, length(considered_resolutions), Q,
                   amp, nrep, overlap_pct, paste(gammas, collapse = "_"), sparsity, alpha, partial_u, chromosome_to_select, tune_gamma_indicator)
write_csv(data.frame(all_res_across_resolutions), out_dir)


out_dir <- sprintf("%s/UKB_kelp_snr_M_%d_Q_%d_amp_%.1f_nrep_%d_ol_%.1f_g_%s_sp%.4f_al%.4f_u%d_c%d_tg%d.csv", save_dir, length(considered_resolutions), Q,
                   amp, nrep, overlap_pct, paste(gammas, collapse = "_"), sparsity, alpha, partial_u, chromosome_to_select, tune_gamma_indicator)
write_csv(data.frame(my_snr), out_dir)


