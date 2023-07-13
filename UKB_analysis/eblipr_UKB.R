#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# IMPORTANT: The data is not available; we have applied for it from a source (the UKBiobank) and interested parties can also apply to the same source

# only need these in this script
library(egg)
library(ggplotify)
library(cowplot)
library(grid)

population <- as.character(args[1])
if(is.na(population)) population <- "british_unrelated"

outcomes <- as.character(args[2])
outcomes <- eval(parse(text = outcomes)) 
if(is.na(outcomes)[1]) outcomes <- "height"

alpha <- as.numeric(args[3])
if(is.na(alpha)) alpha <- 0.1

gamma <- as.character(args[4])
gamma <- eval(parse(text = gamma)) 
if(is.na(gamma)) gamma <- alpha / 8

use_evalues <- as.character(args[5])
use_evalues <- eval(parse(text = use_evalues)) 
if(is.na(use_evalues)) use_evalues <- FALSE

use_multiplier <- as.character(args[6])
use_multiplier <- eval(parse(text = use_multiplier)) 
if(is.na(use_multiplier)) use_multiplier <- FALSE

match_matteo_group <- as.character(args[7])
match_matteo_group <- eval(parse(text = match_matteo_group)) 
if(is.na(match_matteo_group)) match_matteo_group <- FALSE

multiplier <- as.character(args[8])
multiplier <- eval(parse(text = multiplier)) 
if(is.na(multiplier)) multiplier <- NA

print("This is the population used in this script:")
print(population)

print("This is the outcome used in this script:")
print(outcomes)

print("This is alpha used in this script:")
print(alpha)


print("This is gamma used in this script:")
print(gamma)

print("Are e-values used in this script?")
print(use_evalues)

print("Are we matching Matteo's Groups?")
print(match_matteo_group)


print("Are mutlipliers (not uniform, not original e-values used in this script:")
print(use_multiplier)

if(use_multiplier) {
  print("This is the multiplier")
  print(multiplier)
}


resolution_labels <- c("single-SNP", "3 kb", "20 kb", "41 kb", "81 kb", "208 kb", "425 kb")

######## Define Defaults, Set Working Directory and Read in Functions ###########

set.seed(2022)


# Set paths 
if(dir.exists("")){
  mydir <- ""
  setwd(mydir)
  
  my_out_dir <- "" 
}


# load all functions
source(paste0(mydir, "utils/load_aux.R"))

outcome_string <- paste0(outcomes, collapse = "_")

# chromosomes to consider
chr <- seq(1, 22)

# resolutions considered
all_resolutions <- seq(0, 6)

# all groups to be considered
cand_groups <- list()
# all snps considered
snps <- c()

# rejected by global null
rej_global_res <- c()

# rejected by knockoff
rej_knockoff <- c()

snp_group_chr_matrix <- c()

min_non_inf_gamma <- c()

######## Read in Data ############

# this will always be run from scratch 
# read in data from Matteo's results
for(r in all_resolutions) {
  
  
  res <- read_process_UKB_lasso(resolution = r, 
                                population = population,
                                alpha = alpha, 
                                gamma = gamma,
                                mydir = mydir, 
                                dirpath_group = "",
                                dirpath_lasso = "",
                                outdir_path = paste0(my_out_dir, "snp_group_w/"), 
                                outcomes = outcomes, 
                                chr = chr, 
                                smallest_possible_gamma = FALSE, 
                                use_multiplier = use_multiplier, 
                                multiplier = multiplier, 
                                match_matteo_group = match_matteo_group)
  
  cand_groups <- c(cand_groups, res$cand_groups_resolution)
  snps <- c(snps, res$snps)
  rej_global_res <- rbind(rej_global_res, data.frame(res$rejected_matrix_ebh_global_res))
  snp_group_chr_matrix <- rbind(snp_group_chr_matrix, res$snp_group_chr_matrix_evals)
  rej_knockoff <- rbind(rej_knockoff, res$rejected_knockoff)
  min_non_inf_gamma <- rbind(min_non_inf_gamma, data.frame(min_gamma = res$min_non_inf_gamma) %>% mutate(Resolution = r))
  
}

saveRDS(cand_groups, file=paste0(my_out_dir, "RDS/", outcome_string,"_", population ,"_cand_groups_all_res_gamma", gamma))
saveRDS(snp_group_chr_matrix, file=paste0(my_out_dir, "RDS/", outcome_string, "_", population ,"_snp_group_chr_matrix_gamma", gamma))
saveRDS(snps, file=paste0(my_out_dir, "RDS/", outcome_string, "_", population ,"_snps_gamma", gamma))
saveRDS(rej_global_res, file=paste0(my_out_dir, "RDS/", outcome_string, "_", population , "_rej_global_res_gamma_", gamma))

if(!is.null(dim(rej_knockoff)[1])) {
  colnames(rej_knockoff) <- c("CHR", "Group","Resolution","BP.min" , "BP.max","rejected_knockoff")
  saveRDS(rej_knockoff, file=paste0(my_out_dir, "RDS/", outcome_string, "_rej_knockoff_res_gamma_", gamma))
}

########## Basic Data Wrangling ##############

# create chr-group-resolution identifier
snp_group_chr_matrix$Res_CHR_Group <- paste0(snp_group_chr_matrix$Resolution, "_", snp_group_chr_matrix$CHR_Group)

#  group size for each resolution-CHR-group
counts_per_res_chr_group <- snp_group_chr_matrix %>% 
  group_by(Res_CHR_Group) %>%
  count()

# merge group size into dataset
snp_group_chr_matrix_wcounts <- dplyr::left_join(snp_group_chr_matrix, 
                                                 counts_per_res_chr_group, by = c("Res_CHR_Group"))

# calculate weight-adjusted fraction
snp_group_chr_matrix_wcounts <- snp_group_chr_matrix_wcounts %>% 
  mutate(frac_adjusted = fracs / n)

# filter fractions > 0 
snp_group_chr_matrix_wcounts_filtered_nonzero <- snp_group_chr_matrix_wcounts %>% filter(fracs > 0)

if(!is.null(dim(rej_knockoff)[1]) & length(outcomes) > 1 & length(rej_knockoff) > 0) {
  data.frame(rej_knockoff) %>% 
    filter(rejected_knockoff == 1) %>%
    mutate(Res_CHR_Group = paste0(Resolution, "_", CHR, "_", Group)) %>%
    pull(Res_CHR_Group) %>% unique() -> knockoff_rejected_groups
}

# get total number of snps in the groups with non-zero fractions
snps_filtered_eblipr <- snp_group_chr_matrix_wcounts_filtered_nonzero %>% 
  pull(SNP) %>% unique()

# get total number of groups(original cand_groups)
total_groups <- length(cand_groups)

# calculate number of snps in each group for groupsize
number_snp_per_group <- snp_group_chr_matrix %>% 
  group_by(Res_CHR_Group) %>% 
  summarise(number_unique_SNP_per_group = n_distinct(SNP))

saveRDS(number_snp_per_group, file=paste0(my_out_dir, "RDS/", outcome_string, "_", population , "_group_size"))

# calculate average number of SNPs per resolution 
number_snp_per_group_and_res <- snp_group_chr_matrix %>% 
  group_by(Resolution, Res_CHR_Group) %>% 
  summarise(number_unique_SNP_per_group = n_distinct(SNP))

avg_number_snp_per_group_res <- number_snp_per_group_and_res %>% 
  group_by(Resolution) %>% 
  summarise(avg_number_snps_per_group = mean(number_unique_SNP_per_group))

print(avg_number_snp_per_group_res)


########## Plot Global Results ###########

# note that this does not match the knockoff output even for a single outcome 
# because the global results are calculated using the e-value, which is based on gamma 

if(length(outcomes) > 1) {
  
  # get fractions and other info for rejected groups 
  rej_global_groups <- rej_global_res %>% mutate(Res_CHR_Group = paste0(Resolution, "_",  CHR, "_", Group)) %>% 
    pull(Res_CHR_Group) %>% unique()
  rej_global_res_with_frac <- snp_group_chr_matrix %>% filter(Res_CHR_Group %in% rej_global_groups)
  
  
  saveRDS(rej_global_res_with_frac, file=paste0(my_out_dir, "RDS/", outcome_string, "_rej_global_gamma_", gamma))
  
  # plot results separately for each chromosome
  for(c in chr) {
    
    rej_chromosome <- rej_global_res_with_frac %>% filter(CHR == c)
    global_null_plot_chr <- plot_chicago(window.chr = c, 
                                         window.left = min(rej_chromosome$BP.min, na.rm = TRUE),
                                         window.right = max(rej_chromosome$BP.max, na.rm = TRUE),
                                         Discoveries = rej_chromosome, 
                                         evals = FALSE)
    
    ggsave(paste0(my_out_dir, "figures/", outcome_string, "_global_rej_chr_", c, "_gamma_", gamma, ".png"), 
           plot = global_null_plot_chr, device = png, width = 20, height = 5)
    
  }
  
  
  # total number of global rejections by resolution
  rej_global_res_with_frac %>% 
    dplyr::select(Resolution, Res_CHR_Group) %>% 
    unique() %>% 
    group_by(Resolution) %>%
    count() -> total_number_global_rejected_by_resolution
  
  group_res_global <- rej_global_res_with_frac %>% 
    dplyr::select(Res_CHR_Group, Resolution) %>%
    unique() 
  
  inverse_size_rejections_global <- number_snp_per_group %>% 
    filter(Res_CHR_Group %in% unique(rej_global_res_with_frac$Res_CHR_Group)) %>%
    mutate(inverse_size = 1 / number_unique_SNP_per_group) 
  
  group_res_global <- left_join(group_res_global, inverse_size_rejections_global, by = c("Res_CHR_Group"))
  
  power_global_by_resolution <- group_res_global %>% 
    group_by(Resolution) %>% 
    summarise(n = sum(inverse_size)) %>% 
    mutate(n = round(n, 2))
  
  number_global_snps_implicated_by_resolution <- group_res_global %>% 
    group_by(Resolution) %>% 
    summarise(n = sum(number_unique_SNP_per_group)) 
  
  # make sure all resolutions are present 
  number_global_snps_implicated_by_resolution <- left_join(data.frame(Resolution = all_resolutions), number_global_snps_implicated_by_resolution, by = c("Resolution"))
  
  print(power_global_by_resolution)
  
  
  # also filter to outer nodes
  # filter to outer nodes
  outer_nodes_global <- filter_outer_node(rej_global_res_with_frac)
  
  # count the number of rejected outer nodes per resolution
  outer_nodes_global %>% 
    group_by(Resolution) %>% 
    count() -> number_global_rejections_by_resolution 
  
  # plot the rejected outer nodes only
  more_info_outer_nodes_global <- rej_global_res_with_frac %>% filter(Res_CHR_Group %in% outer_nodes_global$Res_CHR_Group)
  
  saveRDS(more_info_outer_nodes_global, file=paste0(my_out_dir, "RDS/", outcome_string, "_rej_outer_global_gamma_", gamma))
  
  # calculate resolution-adjusted "power" 
  inverse_size_rejections_outer_global <- number_snp_per_group %>% 
    filter(Res_CHR_Group %in% unique(more_info_outer_nodes_global$Res_CHR_Group)) %>%
    mutate(inverse_size = 1 / number_unique_SNP_per_group) 
  
  power_global_outer <- sum(inverse_size_rejections_outer_global$inverse_size)
  
  number_global_outer_snps_implicated_by_res <- number_snp_per_group_and_res %>% 
    filter(Res_CHR_Group %in% unique(more_info_outer_nodes_global$Res_CHR_Group)) %>%
    group_by(Resolution) %>% 
    summarise(number_global_outer_snps_implicated = sum(number_unique_SNP_per_group)) 
  
  # make sure all resolutions are present
  number_global_outer_snps_implicated_by_res <- left_join(data.frame(Resolution = all_resolutions),
                                                          number_global_outer_snps_implicated_by_res, 
                                                          by = c("Resolution"))
   
  
  
  
  
  # power by resolution
  
  for(c in chr) {
    
    rej_chromosome <- more_info_outer_nodes_global %>% filter(CHR == c)
    knockoff_rej_plot_chr <- plot_chicago(window.chr = c, 
                                          window.left = min(rej_chromosome$BP.min, na.rm = TRUE),
                                          window.right = max(rej_chromosome$BP.max, na.rm = TRUE),
                                          Discoveries = rej_chromosome, 
                                          evals = FALSE)
    
    ggsave(paste0(my_out_dir, "figures/", outcome_string, "_global_outer_rej_chr_", c, "_gamma_", gamma, ".png"),
           plot = knockoff_rej_plot_chr, device = png, width = 20, height = 5)
    
  }
  
  
  
  
  
  
}

########## Plot Knockoff Rejections  ###########

if(length(outcomes) == 1) {
  
  # get more info for groups rejected by the knockoff
  rej_knockoff_groups <- data.frame(rej_knockoff) %>% 
    filter(rejected_knockoff == 1) %>% 
    mutate(Res_CHR_Group = paste0(Resolution, "_", CHR, "_", Group)) %>%  pull(Res_CHR_Group) %>% unique()
  rej_knockoff_res <- snp_group_chr_matrix %>% filter(Res_CHR_Group %in% rej_knockoff_groups)
  
  saveRDS(rej_knockoff_res, file=paste0(my_out_dir, "RDS/", outcome_string, "_rej_knockoff_info_gamma_", gamma))
  
  # plot results
  for(c in chr) {
    
    rej_chromosome <- rej_knockoff_res %>% filter(CHR == c)
    knockoff_rej_plot_chr <- plot_chicago(window.chr = c, 
                                          window.left = min(rej_chromosome$BP.min, na.rm = TRUE),
                                          window.right = max(rej_chromosome$BP.max, na.rm = TRUE),
                                          Discoveries = rej_chromosome, 
                                          evals = FALSE)
    
    ggsave(paste0(my_out_dir, "figures/", outcome_string, "_knockoff_rej_chr_", c, "_gamma_", gamma, ".png"),
           plot = knockoff_rej_plot_chr, device = png, width = 20, height = 5)
    
  }
  
  # filter to outer nodes
  outer_nodes_knockoff <- filter_outer_node(rej_knockoff_res)
  
  
  # plot the rejected outer nodes only
  more_info_outer_nodes <- rej_knockoff_res %>% filter(Res_CHR_Group %in% outer_nodes_knockoff$Res_CHR_Group)
  
  saveRDS(more_info_outer_nodes, file=paste0(my_out_dir, "RDS/", outcome_string, "_rej_knockoff_outer_gamma_", gamma))
  
  
  for(c in chr) {
    
    rej_chromosome <- more_info_outer_nodes %>% filter(CHR == c)
    knockoff_rej_plot_chr <- plot_chicago(window.chr = c, 
                                          window.left = min(rej_chromosome$BP.min, na.rm = TRUE),
                                          window.right = max(rej_chromosome$BP.max, na.rm = TRUE),
                                          Discoveries = rej_chromosome, 
                                          evals = FALSE)
    
    ggsave(paste0(my_out_dir, "figures/", outcome_string, "_knockoff_outer_rej_chr_", c, "_gamma_", gamma, ".png"),
           plot = knockoff_rej_plot_chr, device = png, width = 20, height = 5)
    
  }
  
  # count the number of rejected outer nodes per resolution
  outer_nodes_knockoff %>% 
    group_by(Resolution) %>% 
    count() -> number_knockoff_outer_rejections_by_resolution 
  
  # count the total number of rejections per resolution by the knockoff 
  rej_knockoff_res %>% 
    dplyr::select(Resolution, Res_CHR_Group) %>% 
    unique() %>% 
    group_by(Resolution) %>% 
    count() -> total_number_knockoff_rejections_by_resolution
  
  # calculate resolution-adjusted "power" 
  inverse_size_rejections_knockoff <- number_snp_per_group %>% 
    filter(Res_CHR_Group %in% unique(outer_nodes_knockoff$Res_CHR_Group)) %>%
    mutate(inverse_size = 1 / number_unique_SNP_per_group) 
  
  power_knockoff_outer <- sum(inverse_size_rejections_knockoff$inverse_size)
  
  number_knockoff_outer_snps_implicated_by_res <- number_snp_per_group_and_res %>% 
    filter(Res_CHR_Group %in% unique(outer_nodes_knockoff$Res_CHR_Group)) %>%
    group_by(Resolution) %>% 
    summarise(number_knockoff_outer_snps_implicated = sum(number_unique_SNP_per_group))
  
  
  # calculate resolution-adjusted "power" for regular knockoff at each level of resolution
  
  inverse_size_rejections_regular_knockoff <- number_snp_per_group %>% 
    filter(Res_CHR_Group %in% unique(rej_knockoff_groups)) %>%
    mutate(inverse_size = 1 / number_unique_SNP_per_group) 
  
  group_res_knockoff <- data.frame(rej_knockoff) %>% 
    filter(rejected_knockoff == 1) %>% 
    mutate(Res_CHR_Group = paste0(Resolution, "_", CHR, "_", Group)) %>%
    dplyr::select(Res_CHR_Group, Resolution) %>%
    unique() 
  
  group_res_knockoff <- left_join(group_res_knockoff, inverse_size_rejections_regular_knockoff, by = c("Res_CHR_Group"))
  
  power_knockoff_by_resolution <- group_res_knockoff %>% 
    group_by(Resolution) %>% 
    summarise(n = sum(inverse_size))
  
  snps_implicated_knockoff_by_resolution <-  group_res_knockoff %>% 
    group_by(Resolution) %>% 
    summarise(n = sum(number_unique_SNP_per_group))
  
  print(power_knockoff_by_resolution)
  
}


############ KELP #############

cand_groups_for_eblipr <- cand_groups[sapply(cand_groups, function(x) x$fracs > 0)]

if(length(cand_groups_for_eblipr) > 0) {
  
  cand_groups_for_eblipr <- lapply(cand_groups_for_eblipr,
                                   function(x) c(x, Res_CHR_Group =paste0(x$Resolution, "_", x$CHR_Group)))
  
  # number groups in optimization procedure
  n_groups_filtered <- length(cand_groups_for_eblipr)
  
  # number of unique snps 
  n_individual_hypotheses_filtered <- length(snps_filtered_eblipr)
  
  # run optimization problem
  old_time <- Sys.time() # get start time
  
  res <- kelp(group_candidates_list = cand_groups_for_eblipr,
              individual_snps_candidates = snps_filtered_eblipr,
              single_u = TRUE,
              alpha = alpha, 
              verbose = T,
              weighted = T, 
              fractions_only = !use_evalues, 
              M = 7, # number of resolutions
              n_individual_hypotheses_filtered = n_individual_hypotheses_filtered , 
              n_groups_filtered = n_groups_filtered, 
              total_groups = total_groups, 
              Q = length(outcomes), 
              u_seq = 1)
  new_time <- Sys.time() - old_time # calculate difference
  print("RUNTIME EBLIPR:")
  print(new_time) 
  eblipr_detections <- res[[1]]
  
  saveRDS(eblipr_detections, paste0(my_out_dir, "RDS/", outcome_string, "_", population, "eblipr_detections_gamma", gamma))
  
  # get selected groups
  selected_eblipr_groups <- sapply(eblipr_detections, `[[`, "Res_CHR_Group")
  
  print("Number of rejections by eblipr")
  print(length(selected_eblipr_groups))
  
  # Check that all SNPs are unique
  selected_eblipr_groups_SNPs <- unlist(sapply(eblipr_detections, `[[`, "SNP"))
  length(unique(selected_eblipr_groups_SNPs)) == length(selected_eblipr_groups_SNPs)
  
  rejected_eblipr_with_info <- snp_group_chr_matrix_wcounts_filtered_nonzero %>% 
    dplyr::filter(Res_CHR_Group %in% selected_eblipr_groups)
  
  ####### Plot Eblipr Results ####
  
  # to generate plot find all groups in lower resolutions implicitly included
  additional_groups_implicitly_rejected <- c()
  
  # need to get all groups in higher-up resolutions 
  for(g in selected_eblipr_groups) {
    
    # resolution rejected group is in
    resolution_of_group_g <- snp_group_chr_matrix %>% dplyr::filter(Res_CHR_Group == g) %>% pull(Resolution) %>% unique()
    
    # snps that are included in the rejected group
    SNPs_in_group_g <- snp_group_chr_matrix %>% dplyr::filter(Res_CHR_Group == g) %>% pull(SNP) %>% unique()
    
    # find all groups in lower resolutions (higher resolution number) that contain these snps
    higher_resolution_groups_containing_g <- snp_group_chr_matrix %>% 
      dplyr::filter(SNP %in% SNPs_in_group_g) %>% 
      dplyr::filter(Resolution > resolution_of_group_g) %>%
      pull(Res_CHR_Group) %>% unique()
    
    additional_groups_implicitly_rejected <- c(additional_groups_implicitly_rejected, higher_resolution_groups_containing_g)
    
  }
  
  # save results
  saveRDS(additional_groups_implicitly_rejected, 
          paste0(my_out_dir, "RDS/", outcome_string, "_", population, "eblipr_add_groups", gamma))
  saveRDS(rejected_eblipr_with_info, 
          paste0(my_out_dir, "RDS/", outcome_string, "_", population, "eblipr_detections_info_gamma", gamma))
  
  more_info_additional_groups_implicity_rejected <- snp_group_chr_matrix %>% filter(Res_CHR_Group %in% additional_groups_implicitly_rejected)
  saveRDS(more_info_additional_groups_implicity_rejected, paste0(my_out_dir, "RDS/", outcome_string, "_", population, 
                                                                 "eblipr_add_info_gamma", gamma))
  
  
  
  for(c in chr) {
    
    more_info_chromosome_additional_ones <- more_info_additional_groups_implicity_rejected %>% filter(CHR == c)
    
    rej_chromosome <- rejected_eblipr_with_info %>% filter(CHR == c)
    eblipr_plot_chr <- plot_chicago_with_bigger_res(window.chr = c, 
                                                    window.left = min(rej_chromosome$BP.min, na.rm = TRUE),
                                                    window.right = max(rej_chromosome$BP.max, na.rm = TRUE),
                                                    Discoveries = rej_chromosome, 
                                                    more_info_additional_groups_implicity_rejected = more_info_chromosome_additional_ones,
                                                    evals = FALSE)
    
    ggsave(paste0(my_out_dir, "figures/", outcome_string, "_eblipr_chr", c, "_gamma_" , gamma,".png"),
           plot = eblipr_plot_chr, device = png, width = 20, height = 5)
    
    
    eblipr_plot_single_line_chr <- plot_chicago_single_line(window.chr = c, 
                                                            window.left = min(rej_chromosome$BP.min, na.rm = TRUE),
                                                            window.right = max(rej_chromosome$BP.max, na.rm = TRUE),
                                                            Discoveries = rej_chromosome, 
                                                            evals = FALSE)
    
    length_plot <- round((max(rej_chromosome$BP.max, na.rm = TRUE) - min(rej_chromosome$BP.min, na.rm = TRUE))*1e-6, 3)/5
    
    ggsave(paste0(my_out_dir, "figures/", outcome_string, "_eblipr_chr_single_line", c, "_gamma_" , gamma,".png"),
           plot = eblipr_plot_single_line_chr, device = png, width = length_plot, height = 5)
    
    
    
    
    functional_annotation_with_barplot <- plot_combined(window.chr = c, 
                                                        window.left = min(rej_chromosome$BP.min, na.rm = TRUE),
                                                        window.right = max(rej_chromosome$BP.max, na.rm = TRUE),
                                                        Discoveries = rej_chromosome, 
                                                        data_dir = paste0(mydir, "UKB_data_input"),
                                                        highlight.gene=NULL, max.gene.rows=10)
    
    
    
    ggsave(paste0(my_out_dir, "figures/", outcome_string, "_eblipr_chr_fa_single_line_chr", c, "_gamma_" , gamma,".png"),
           plot = functional_annotation_with_barplot, device = png, width = 20, height = 5)
  }
  
  
  
  # save number of rejections for each resolution by knockoff and by eblipr
  rejected_eblipr_with_info %>% 
    dplyr::select(Resolution, Res_CHR_Group) %>% 
    unique() %>% 
    group_by(Resolution) %>%
    count() -> number_eblipr_rejected_by_resolution
  
  
  # calculate resolution-adjusted "power" 
  inverse_size_rejections_eblipr <- number_snp_per_group %>% 
    filter(Res_CHR_Group %in% unique(selected_eblipr_groups)) %>%
    mutate(inverse_size = 1 / number_unique_SNP_per_group) 
  
  number_eblipr_snps_implicated_by_res <- number_snp_per_group_and_res %>% 
    filter(Res_CHR_Group %in% unique(selected_eblipr_groups)) %>%
    group_by(Resolution) %>% 
    summarise(number_eblipr_snps_implicated = sum(number_unique_SNP_per_group))
  
  power_eblipr <- sum(inverse_size_rejections_eblipr$inverse_size)
  
  ############### focused e-BH ################# 
  
  cat("running focusedeBH ...")
  group_fracs <- snp_group_chr_matrix %>% dplyr::select(Res_CHR_Group, 
                                                        fracs, 
                                                        Resolution,
                                                        CHR_Group, 
                                                        Group, 
                                                        CHR) %>% unique()
  
  rejections_focusedeBH <- focusedeBH_UKB(unique_groups = group_fracs, 
                                          snp_groups = snp_group_chr_matrix,
                                          numerator_threshold = 7, 
                                          alpha = 0.1) 
  rejections_focusedeBH %>%
    group_by(Resolution) %>%
    count() -> rejections_focusedeBH_by_resolution
  
  
  # calculate resolution-adjusted "power" 
  inverse_size_rejections_fbh <- number_snp_per_group %>% 
    filter(Res_CHR_Group %in% rejections_focusedeBH$Res_CHR_Group) %>%
    mutate(inverse_size = 1 / number_unique_SNP_per_group) 
  
  number_focusedeBH_snps_implicated_by_res <- number_snp_per_group_and_res %>% 
    filter(Res_CHR_Group %in% rejections_focusedeBH$Res_CHR_Group) %>%
    group_by(Resolution) %>% 
    summarise(number_focusedeBH_snps_implicated = sum(number_unique_SNP_per_group))
  
  power_focusedeBH <- sum(inverse_size_rejections_fbh$inverse_size)
  
  cat("done.\n")
  saveRDS(rejections_focusedeBH, file=paste0(my_out_dir, "RDS/", outcome_string, "_", population , "_rej_fbh_gamma_", gamma))
  print("rejections febh:")
  print(dim(rejections_focusedeBH))
  
  
  
  # compare focused eBH rejections to eblipr 
  fbh_not_eblipr <- setdiff(unique(rejections_focusedeBH$Res_CHR_Group), 
                            unique(rejected_eblipr_with_info$Res_CHR_Group))
  
  print("fbh not eblipr")
  print(fbh_not_eblipr)
  eblipr_not_fbh <- setdiff(unique(rejected_eblipr_with_info$Res_CHR_Group), 
                            unique(rejections_focusedeBH$Res_CHR_Group))
  print("eblipr not fbh")
  print(eblipr_not_fbh)
  
  if(length(fbh_not_eblipr) > 0 | length(eblipr_not_fbh) > 0) {
    fbh_not_eblipr <- data.frame(Res_CHR_Group = fbh_not_eblipr, Method = "fbh_not_eblipr") 
    fbh_not_eblipr <- left_join(fbh_not_eblipr, number_snp_per_group, by = c("Res_CHR_Group"))
    
    
    
    eblipr_not_fbh <- data.frame(Res_CHR_Group = eblipr_not_fbh, Method = "eblipr_not_fbh")
    eblipr_not_fbh <- left_join(eblipr_not_fbh, number_snp_per_group, by = c("Res_CHR_Group"))
    
    
    not_matching_fbh_eblipr_groups <- rbind(fbh_not_eblipr,eblipr_not_fbh )
    saveRDS(not_matching_fbh_eblipr_groups, file=paste0(my_out_dir, "RDS/", outcome_string, "_", population , "_not_match_fbh_eblipr_groups_", gamma))
    
    fbh_not_eblipr_group_snps <- snp_group_chr_matrix %>% filter(Res_CHR_Group %in% c(fbh_not_eblipr$Res_CHR_Group))
    eblipr_not_fbh_group_snps <- snp_group_chr_matrix %>% filter(Res_CHR_Group %in% c(eblipr_not_fbh$Res_CHR_Group))
    
    snps_fbh_not_eblipr <- setdiff(unique(fbh_not_eblipr_group_snps$SNP), unique(eblipr_not_fbh_group_snps$SNP))
    snps_eblipr_not_fbh <- setdiff(unique(eblipr_not_fbh_group_snps$SNP), unique(fbh_not_eblipr_group_snps$SNP))
    
    snps_not_matching_fbh_eblipr <- c(snps_eblipr_not_fbh, snps_fbh_not_eblipr)
    saveRDS(snps_not_matching_fbh_eblipr, file=paste0(my_out_dir, "RDS/", outcome_string, "_", population , "_snps_not_match_fbh_eblipr_", gamma))
    
    saveRDS(fbh_not_eblipr_group_snps, file=paste0(my_out_dir, "RDS/", outcome_string, "_", population , "_fbh_not_match_eblipr_group_snps_", gamma))
    saveRDS(eblipr_not_fbh_group_snps, file=paste0(my_out_dir, "RDS/", outcome_string, "_", population , "_eblipr_not_match_fbh_group_snps_", gamma))
    
    
    eblipr_not_fbh_group_snps <- eblipr_not_fbh_group_snps %>% 
      mutate(Res_CHR_Group_kelp_not_fbh = Res_CHR_Group) %>% 
      dplyr::select(SNP,Res_CHR_Group_kelp_not_fbh) 
    
    
    fbh_not_eblipr_group_snps <- fbh_not_eblipr_group_snps %>% 
      mutate(Res_CHR_Group_fbh_not_kelp = Res_CHR_Group) %>%
      dplyr::select(SNP,Res_CHR_Group_fbh_not_kelp) 
    
    
    merged_snp_group <- inner_join(fbh_not_eblipr_group_snps, eblipr_not_fbh_group_snps, by = c("SNP"))
    merged_snp_group <- merged_snp_group %>% dplyr::select(Res_CHR_Group_fbh_not_kelp, Res_CHR_Group_kelp_not_fbh) %>% unique()
    
    saveRDS(merged_snp_group, file=paste0(my_out_dir, "RDS/", outcome_string, "_", population , "_same_groups_not_match_", gamma))
    
  }
  
  
  # combine into one data frame and save
  if(length(outcomes) == 1) {
    number_rejections_by_resolution_df <- rbind(number_eblipr_rejected_by_resolution %>% mutate(Method = "Kelp"), 
                                                rejections_focusedeBH_by_resolution %>% mutate(Method = "Focused eBH"),
                                                number_knockoff_outer_rejections_by_resolution  %>% mutate(Method = "Outer Nodes"), 
                                                total_number_knockoff_rejections_by_resolution  %>% mutate(Method = "All"), 
                                                power_knockoff_by_resolution %>% mutate(Method = "Power"), 
                                                snps_implicated_knockoff_by_resolution %>% mutate(Method = "SNPs"))
    
    number_rejections_by_resolution_df <- number_rejections_by_resolution_df %>% dplyr::select(Method, Resolution, n) 
    number_rejections_by_resolution_df_wide <- number_rejections_by_resolution_df %>% 
      pivot_wider(names_from = Method, values_from = n) %>% 
      mutate(Difference = `Outer Nodes` - Kelp)
    
    
    print(sum(number_rejections_by_resolution_df_wide$Difference) / sum(number_rejections_by_resolution_df_wide$`Outer Nodes`))
    
    
    number_rejections_by_resolution_df_wide <- data.frame(number_rejections_by_resolution_df_wide) 
    
    powers <- data.frame(power = c(power_focusedeBH, power_knockoff_outer, power_eblipr), 
                         method = c("focused eBH", "outer", "kelp"))
    
    snps_implicated_by_res <- cbind(number_eblipr_snps_implicated_by_res$Resolution,
                                    number_eblipr_snps_implicated_by_res$number_eblipr_snps_implicated,
                                    number_knockoff_outer_snps_implicated_by_res$number_knockoff_outer_snps_implicated, 
                                    snps_implicated_knockoff_by_resolution$n)
    total_snps_implicated <- c("Total", sum(number_eblipr_snps_implicated_by_res$number_eblipr_snps_implicated), 
                               sum(number_knockoff_outer_snps_implicated_by_res$number_knockoff_outer_snps_implicated), 
                               NA)
    
    colnames(snps_implicated_by_res) <- c("Resolution", "KeLP_snps", "outer_snps", "all_snps")
    
  } else {
    number_rejections_by_resolution_df <- rbind(number_eblipr_rejected_by_resolution %>% mutate(Method = "Kelp"), 
                                                rejections_focusedeBH_by_resolution %>% mutate(Method = "Focused eBH"),
                                                number_global_rejections_by_resolution %>% mutate(Method = "Outer Nodes"), 
                                                total_number_global_rejected_by_resolution %>% mutate(Method = "All"), 
                                                power_global_by_resolution %>% mutate(Method = "Power"), 
                                                number_global_snps_implicated_by_resolution %>% mutate(Method = "SNPs"))
    
    number_rejections_by_resolution_df <- number_rejections_by_resolution_df %>% dplyr::select(Method, Resolution, n)
    
    number_rejections_by_resolution_df <- number_rejections_by_resolution_df %>% dplyr::select(Method, Resolution, n) 
    
    number_rejections_by_resolution_df_wide <- number_rejections_by_resolution_df %>% 
      pivot_wider(names_from = Method, values_from = n) %>% 
      mutate(Difference = `Outer Nodes` - Kelp)
    
    print(sum(number_rejections_by_resolution_df_wide$Difference) / sum(number_rejections_by_resolution_df_wide$`Outer Nodes`))
    
    
    powers <- data.frame(power = c(power_focusedeBH, power_global_outer, power_eblipr), 
                         method = c("focused eBH", "outer", "kelp"))
    
    
    snps_implicated_by_res <- cbind(number_eblipr_snps_implicated_by_res$Resolution,
                                    number_eblipr_snps_implicated_by_res$number_eblipr_snps_implicated,
                                    number_global_outer_snps_implicated_by_res$number_global_outer_snps_implicated, 
                                    number_global_snps_implicated_by_resolution$n)
    total_snps_implicated <- c("Total", sum(number_eblipr_snps_implicated_by_res$number_eblipr_snps_implicated), 
                               sum(number_global_outer_snps_implicated_by_res$number_global_outer_snps_implicated, na.rm = TRUE), 
                               NA)
    colnames(snps_implicated_by_res) <- c("Resolution", "KeLP_snps", "outer_snps", "all_snps")
    
  }
  
  # save as RDS 
  saveRDS(number_rejections_by_resolution_df_wide, file=paste0(my_out_dir, "RDS/", outcome_string, "_", population , "_rej_by_res_table_", gamma))
  
  print(number_rejections_by_resolution_df_wide)
  
  if(length(outcomes) == 1) {
    number_rejections_by_resolution_df_wide <- number_rejections_by_resolution_df_wide %>% 
      dplyr::select(Resolution, Kelp, Outer.Nodes, All)
    colnames(number_rejections_by_resolution_df_wide) <- c("Resolution", "KeLP", "outer nodes", "All")
    
    # add total! 
    
    total_number_of_rejections <- data.frame(Resolution = "Total", 
                                             KeLP = sum(number_rejections_by_resolution_df_wide$KeLP, na.rm = TRUE), 
                                             `outer nodes` = sum(number_rejections_by_resolution_df_wide$`outer nodes`, na.rm = TRUE), 
                                             All = "")
    

    colnames(total_number_of_rejections) <- c("Resolution", "KeLP", "outer nodes", "All")
  } else {
    number_rejections_by_resolution_df_wide <- number_rejections_by_resolution_df_wide %>% 
      dplyr::select(Resolution, Kelp, `Outer Nodes`, `All`)
    colnames(number_rejections_by_resolution_df_wide) <- c("Resolution", "KeLP", "outer nodes", "All")
    
    
    # add total! 
    total_number_of_rejections <- data.frame(Resolution = "Total", 
                                             KeLP = sum(number_rejections_by_resolution_df_wide$KeLP, na.rm = TRUE), 
                                             `outer nodes` = sum(number_rejections_by_resolution_df_wide$`outer nodes`, na.rm = TRUE), 
                                             All = NA)

    colnames(total_number_of_rejections) <- c("Resolution", "KeLP", "outer nodes", "All")
    
  }
  
  
  number_rejections_by_resolution_df_wide$Resolution <- resolution_labels
  number_rejections_by_resolution_df_wide[is.na(number_rejections_by_resolution_df_wide)] <- 0
  
  snps_implicated_by_res[is.na(snps_implicated_by_res)] <- 0
  
  number_rejections_by_resolution_df_wide <- rbind(number_rejections_by_resolution_df_wide, total_number_of_rejections)
  
  # combine with SNPs implicated 
  snps_implicated_by_res <- data.frame(snps_implicated_by_res)
  snps_implicated_by_res$Resolution <- resolution_labels
  snps_implicated_by_res <- rbind(snps_implicated_by_res, total_snps_implicated)
  
  
  # merge both together! 
  final_UKB_table <- inner_join(number_rejections_by_resolution_df_wide, snps_implicated_by_res, by = "Resolution")
  
  final_UKB_table$empty_column_1 <- ""
  final_UKB_table$empty_column_2 <- ""
  
  final_UKB_table <- final_UKB_table %>% dplyr::select(Resolution, 
                                                       empty_column_1, 
                                                       KeLP, `outer nodes`, All,
                                                       empty_column_2,
                                                       KeLP_snps, outer_snps, all_snps)
  
  
  colnames(final_UKB_table) <- c("Resolution", "empty","KeLP", "Outer nodes", "Res. spec.","empty","KeLP", "Outer nodes", "Res. spec.")
  
  
  final_UKB_table[,2:9] <- sapply(final_UKB_table[,2:9],as.numeric)
  
  stargazer(final_UKB_table,
            type = "latex",
            summary=FALSE, rownames=FALSE,
            digits = 2, 
            decimal.mark = ",") -> table_res
  
  vars = c("empty.1" = "empty", "KeLP.1" = "KeLP", "Outer nodes.1" = "Outer nodes", "Res. spec..1" = "Res. spec.")
  table_res <- gsub(paste(names(vars), collapse = " & "), paste(vars, collapse = " & "), table_res)
  
  
  tabular_positions <- grep("tabular", table_res)
  table_res <- table_res[tabular_positions[1]:tabular_positions[2]]
  tabular_positions <- grep("tabular", table_res)
  table_res <- table_res[(tabular_positions[1]+1):(tabular_positions[2]-1)]
  table_res <- c(table_res, "\\insertTableNotes")
  table_res[3] <- gsub('\\.', ' ', table_res[3])
  table_res[3] <- gsub('empty', '', table_res[3])
  
  
  part1 <- table_res[1:2]
  
  part2a <- "& & \\multicolumn{3}{c}{Number of rejections} & & \\multicolumn{3}{c}{SNPs implicated}  \\\\ "
  part2b <- "\\cline{3-5} \\cline{7-9} \\\\"
  part3 <- table_res[3:13]
  
  final_table <- c(part1, part2a, part2b, part3)
  
  final_table_p1 <- final_table[1:13]
  final_table_p2 <- final_table[14:15]
  table_res_insert <- "\\hline \\\\[-1.8ex] " 
  table_res <- c(final_table_p1, table_res_insert, final_table_p2)
  table_res <- c(table_res, "\\insertTableNotes")
  table_res <- gsub('NA', ' ', table_res) # ADDED 415
  write(table_res, file=paste0(my_out_dir, "tables/", outcome_string,  "_",population ,"_UKB_comparison.tex"))
  
  
  
  
  
  
  # export power table
  print(powers)
  stargazer(powers,
            type = "latex",
            summary=FALSE, rownames=FALSE,
            digits = 2) -> table_res
  
  tabular_positions <- grep("tabular", table_res)
  table_res <- table_res[tabular_positions[1]:tabular_positions[2]]
  tabular_positions <- grep("tabular", table_res)
  table_res <- table_res[(tabular_positions[1]+1):(tabular_positions[2]-1)]
  table_res <- c(table_res, "\\insertTableNotes")
  table_res[3] <- gsub('\\.', ' ', table_res[3])
  
  write(table_res, file=paste0(my_out_dir, "tables/", outcome_string,  "_",population ,"_power_comparison.tex"))
  

  
}



