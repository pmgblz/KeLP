######################################################################

# Compare the SNPs in Yengo (2022) to the Knockoff Outer Node Rejections #


# This file does the following:
# 1. Read in knockoff outer node and Yengo (2022) rejections
# 2. Find matches for knockoff outer node rejections by name and location 
# 3. Compare distance of not matched groups to next closet SNP

######################################################################


# set paths
if(dir.exists("")){
  mydir <- ""
  setwd(mydir)
  
  
  input_path_yengo <- ""
  my_out_dir <- "" 
}


############ GET SNPS IN YENGO 2022 PAPER ################
# SNPs in Yengo (downloadable from https://www.nature.com/articles/s41586-022-05275-y)
meta_snps_yengo_nature2022 <- read_excel(paste0(input_path_yengo, "supplementary_tables_yengo_nature2022.xlsx"), 
                                         sheet = "ST10 - COJO - METAFE", range = "B2:E12114")


colnames(meta_snps_yengo_nature2022) <- c("CHR", "SNP", "BP.HG19", "BP.HG38")

############ MATCHING THE REJECTIONS ################

# Information on location for snps rejected by eblipr 
height_rejected_knock_with_info <- readRDS(paste0(mydir, "INPUT_DATA_PATH"))

########### ...matching by name #############

# Which snp's can we automatically match by name
height_rejected_knock_with_info <- height_rejected_knock_with_info %>% 
  mutate(direct_match = ifelse(SNP %in% meta_snps_yengo_nature2022$SNP, 1, 0)) 
sum(height_rejected_knock_with_info$direct_match)

# check which groups contain one of the SNPs that could be matched directly
group_contains_direct_match <- height_rejected_knock_with_info %>% 
  group_by(Res_CHR_Group) %>% 
  summarise(group_has_direct_match = sum(direct_match)) %>% 
  mutate(group_has_at_least_one_direct_match = ifelse(group_has_direct_match > 0, 1, 0))
# number of groups containing at least one match by name
sum(group_contains_direct_match$group_has_at_least_one_direct_match)
sum(group_contains_direct_match$group_has_direct_match)

height_rejected_knock_with_info <- left_join(height_rejected_knock_with_info,
                                              group_contains_direct_match, by = "Res_CHR_Group")

# groups that do not contain any snp that can be matched by name
height_rejected_knock_with_info %>% 
  filter(group_has_at_least_one_direct_match == 0) %>% 
  pull(Res_CHR_Group) %>% 
  unique() -> Res_CHR_Group_No_Match

# groups that contain snps that can be matched by name
height_rejected_knock_with_info %>% 
  filter(group_has_at_least_one_direct_match > 0) %>% 
  pull(Res_CHR_Group) %>% 
  unique() -> Res_CHR_Group_Match

# all groups that are rejected 
all_groups <- height_rejected_knock_with_info %>% pull(Res_CHR_Group) %>% unique()

# filter out all groups that cannot be matched by name 
height_rejections_no_name_match <- height_rejected_knock_with_info %>% 
  filter(group_has_at_least_one_direct_match == 0)

# check that there are no inconsistencies 
name_match_snp <- height_rejected_knock_with_info %>% filter(direct_match == 1) %>% dplyr::select(SNP, BP.min, BP.max)
name_match_snp <- left_join(name_match_snp, meta_snps_yengo_nature2022, by = "SNP")


name_match_snp$bigger_than_min <- ifelse(name_match_snp$BP.HG19 >= name_match_snp$BP.min, 1, 0)
name_match_snp$smaller_than_max <- ifelse(name_match_snp$BP.HG19 <= name_match_snp$BP.max, 1, 0)
name_match_snp$in_range <- name_match_snp$bigger_than_min + name_match_snp$smaller_than_max

summary(name_match_snp$in_range)
########### ...matching by location #############

# check if we can we actually match more
corresponding_group <- c()
counter <- 0
for(i in 1:nrow(meta_snps_yengo_nature2022)) {
  counter <- counter + 1
  
  if(counter %% 1000 == 0) {
    print(counter)
  }
  
  snp_row <- meta_snps_yengo_nature2022[i, ]
  
  # find corresponding rejection 
  # only consider groups that do not contain a snp that was previously rejected
  corresponding_rejection <- height_rejections_no_name_match %>% 
    filter(CHR == snp_row$CHR, BP.max >= snp_row$BP.HG19, BP.min <= snp_row$BP.HG19) %>% 
    mutate(yengo_snp = snp_row$SNP, yengo_bp = snp_row$BP.HG19 )
  
  corresponding_group <- rbind(corresponding_group, corresponding_rejection)
  
  
}

# for how many of the rejected groups did we find a match? 
length(unique(corresponding_group$Res_CHR_Group)) 

# for which ones can we not find a match by location and not by name
height_rejections_no_name_no_location_match <- height_rejected_knock_with_info %>% 
  filter(group_has_at_least_one_direct_match == 0) %>% 
  filter(!(Res_CHR_Group %in% corresponding_group$Res_CHR_Group))

length(unique(height_rejections_no_name_no_location_match$Res_CHR_Group))

# how many yengo_snps does this cover? 
length(unique(corresponding_group$yengo_snp))

# which ones remain? 
not_matched_name_location <- height_rejections_no_name_no_location_match %>% 
  dplyr::select(BP.min, BP.max, Res_CHR_Group) %>% unique()
# number of not matched groups
nrow(not_matched_name_location)

# those are groups where no SNP from yengo falls INTO one of the groups  
not_matched_name_location$matched_snp <- ""
not_matched_name_location$matched_bp <- NA
not_matched_name_location$bp_diff <- NA

for(i in 1:nrow(not_matched_name_location)) {
  
  # check distance to min and max BP for each group
  diffs_bp_min <- abs(meta_snps_yengo_nature2022$BP.HG19 - not_matched_name_location[i, ]$BP.min)
  min_diffs_bp_min <- diffs_bp_min[which.min(diffs_bp_min)]
  
  diffs_bp_max <- abs(meta_snps_yengo_nature2022$BP.HG19 - not_matched_name_location[i, ]$BP.min)
  min_diffs_bp_max <- diffs_bp_min[which.min(diffs_bp_max)]
  
  if(min_diffs_bp_min <= min_diffs_bp_max ) {
    min_index <- which.min(diffs_bp_min)
    diffs <- diffs_bp_min
  } else {
    min_index <- which.min(diffs_bp_max)
    diffs <- diffs_bp_max
  }
  
  not_matched_name_location[i, ]$matched_snp <- meta_snps_yengo_nature2022[min_index, ]$SNP
  not_matched_name_location[i, ]$matched_bp <- meta_snps_yengo_nature2022[min_index, ]$BP.HG19
  not_matched_name_location[i, ]$bp_diff <- diffs[min_index]
  
}

# plot histogram or density plot for distance
not_matched_name_location %>% 
  filter(bp_diff <= quantile(bp_diff, 0.95)) %>% 
  ggplot(aes(x = bp_diff)) + geom_histogram(bins = 75) + 
  theme_minimal() + labs(x = "Difference in BP", y = "Count") -> histogram_bottom95_bpdiff_yengo

ggsave(paste0(my_out_dir, "figures/", "histogram_bottom95_bpdiff_yengo_hg19_knockoff.png"), 
       plot = histogram_bottom95_bpdiff_yengo, device = png, width = 7, height = 5)


summary(not_matched_name_location$bp_diff)

# what's the distance between the groups that we CAN match and the closest snp that is not the one that we matched to? 
corresponding_group_snp <- corresponding_group %>% 
  dplyr::select(Res_CHR_Group, yengo_snp) %>% 
  unique() 


# which ones can we match by either name or location
height_rejections_name_location_match <- height_rejected_knock_with_info %>% 
  filter(group_has_at_least_one_direct_match > 0 | Res_CHR_Group %in% corresponding_group$Res_CHR_Group) 

# total number of matched groups
length(unique(height_rejections_name_location_match$Res_CHR_Group))

unique_matching_groups <- unique(height_rejections_name_location_match$Res_CHR_Group)

# find unique snps that we can either match by name or by location
matched_snps <- c()

for(g in unique_matching_groups) {
  
  direct_match <- height_rejections_name_location_match %>%
    filter(Res_CHR_Group == g) %>% 
    pull(group_has_at_least_one_direct_match) %>% 
    unique()
  
  
  if(direct_match == 1) {
    
    # remove all snps in the group, if the group contains a match by name
    matched_snp <- height_rejections_name_location_match %>% filter(Res_CHR_Group == g, direct_match == 1) %>% 
      pull(SNP) %>% 
      unique()
    
  } else {
    
    # get all yengo snps corresponding to the location
    matched_snp <- corresponding_group %>% 
      filter(Res_CHR_Group == g) %>% 
      pull(yengo_snp) %>% 
      unique()
    
  }
  
  
  matched_snps <- c(matched_snps, matched_snp)
  
  
}

# remove all matched snps from the snps we can consider in yengo
# check the distance in location to the next best snp assuming we could not match 
matched_snps <- unique(matched_snps)

# check that the number of snps that are matched makes sense
# total matched snps = snps by location + snps by name
length(matched_snps) == length(unique(corresponding_group$yengo_snp)) + sum(group_contains_direct_match$group_has_direct_match)

matched_name_location <- height_rejections_name_location_match %>% 
  dplyr::select(BP.min, BP.max, Res_CHR_Group) %>% unique()

matched_name_location$next_matched_snp <- ""
matched_name_location$next_matched_bp <- NA
matched_name_location$next_bp_diff <- NA

meta_snps_yengo_nature2022_match_removed <- meta_snps_yengo_nature2022 %>% 
  filter(!(SNP %in% matched_snps))

# how many snps in yengo do our rejections cover
dim(meta_snps_yengo_nature2022)[1] - dim(meta_snps_yengo_nature2022_match_removed)[1]

for(i in 1:nrow(matched_name_location)) {
  
  # check distance to min and max BP for each group
  diffs_bp_min <- abs(meta_snps_yengo_nature2022_match_removed$BP.HG19 - matched_name_location[i, ]$BP.min)
  min_diffs_bp_min <- diffs_bp_min[which.min(diffs_bp_min)]
  
  diffs_bp_max <- abs(meta_snps_yengo_nature2022_match_removed$BP.HG19 - matched_name_location[i, ]$BP.min)
  min_diffs_bp_max <- diffs_bp_min[which.min(diffs_bp_max)]
  
  if(min_diffs_bp_min <= min_diffs_bp_max ) {
    min_index <- which.min(diffs_bp_min)
    diffs <- diffs_bp_min
  } else {
    min_index <- which.min(diffs_bp_max)
    diffs <- diffs_bp_max
  }
  
  matched_name_location[i, ]$next_matched_snp <- meta_snps_yengo_nature2022_match_removed[min_index, ]$SNP
  matched_name_location[i, ]$next_matched_bp <- meta_snps_yengo_nature2022_match_removed[min_index, ]$BP.HG19
  matched_name_location[i, ]$next_bp_diff <- diffs[min_index]
  
}

summary(matched_name_location$next_bp_diff)

# plot histogram or density plot for distance
matched_name_location %>% 
  filter(next_bp_diff <= quantile(next_bp_diff, 0.95)) %>% 
  ggplot(aes(x = next_bp_diff)) + geom_histogram(bins = 75) + 
  theme_minimal() + labs(x = "Difference in BP", y = "Count") -> histogram_bottom95_bpdiff_matched_next_yengo

ggsave(paste0(my_out_dir, "figures/", "histogram_bottom95_bpdiff_matched_next_yengo_hg19_knockoff.png"), 
       plot = histogram_bottom95_bpdiff_matched_next_yengo, device = png, width = 7, height = 5)


# overlap the two histograms
matched_name_location %>% 
  mutate(bp_diff = next_bp_diff) %>%
  dplyr::select(Res_CHR_Group, bp_diff) %>% 
  mutate(Matched = "Yes") %>%
  unique() -> next_snp_for_matched

not_matched_name_location %>% 
  dplyr::select(Res_CHR_Group, bp_diff) %>% 
  mutate(Matched = "No") %>%
  unique() -> next_snp_for_unmatched


next_snps <- rbind(next_snp_for_matched, next_snp_for_unmatched)

next_snps  %>% 
  filter(bp_diff <= quantile(bp_diff, 0.95)) %>% 
  ggplot(aes(x = bp_diff, fill = Matched, colour = Matched)) + 
  geom_density(alpha = 0.5, position = "identity")  + 
  theme_minimal() + labs(x = "Difference in BP", y = "Density") + 
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) -> density_bottom95_bpdiff_next_yengo

ggsave(paste0(my_out_dir, "figures/", "density_bottom95_bpdiff_next_yengo_hg19_knockoff.png"), 
       plot = density_bottom95_bpdiff_next_yengo, device = png, width = 7, height = 5)

next_snps %>% 
  mutate(Method = "Knockoff") -> next_snps 

saveRDS(next_snps, paste0(my_out_dir, "RDS/next_snps_knockoff"))

