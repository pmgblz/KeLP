
###################################################################### 

# Summarize hierarchical simulation results and produce power / fdp plots #

###################################################################### 


library(ggpubr)



if(dir.exists("")){
  mydir <- ""
  setwd(mydir)
  
  my_out_dir <- "" 
}


np_ratio_max = 4

nrep <- 100
alpha <- 0.2    
sparsity = c(0.025, 0.05, 0.1) 
overlap_pct = c(0.5) 
p <- 1000
dpo <- c("0.15") 

#################### ACROSS RESOLUTION VS SEPARAETE BY RESOLUTION ############


# SPARATE PLOTS FOR EACH OVERLAP PCT, WITHIN EACH PLOT ONE COLOR FOR EACH SPARSITY; VARY NP RATIO
# ntg = 1 is gamma / 2
for(ntg in c(1)) {

for(po in dpo) {

for(ol in overlap_pct) {
  
  
  all_res <- data.frame()
  
  for(s in sparsity) {
    
    
    power_fdp <- data.frame()
    
    my_files <- list.files(mydir)
    my_files <- my_files[str_detect(my_files, "hp_union_")]
    my_files <- my_files[!str_detect(my_files, "hp_union_by")]
    my_files <- my_files[str_detect(my_files, paste0("nrep_", nrep, "_"))]
    my_files <- my_files[str_detect(my_files, paste0("sp_", s))]
    my_files <- my_files[str_detect(my_files, paste0("ol_", ol))]
    my_files <- my_files[str_detect(my_files, paste0("dpo", po))]
    my_files <- my_files[str_detect(my_files, paste0("ntg", ntg))]
    my_files <- my_files[str_detect(my_files, paste0("p", p))]
    
    ldf <- lapply(my_files, read.csv)
    
    for(l in 1:length(ldf)) {
      
      
      snr_res_averaged <- ldf[[l]]  %>% 
        dplyr::filter(np_ratio <= np_ratio_max) %>%
        group_by(np_ratio) %>% 
        group_by(method, np_ratio) %>%
        mutate(method = method) %>% 
        mutate(method =  recode_factor(method, eFBH = "e-FBH", knockoffouter = "knockoff outer", ebhknockoffouter = "ebh knockoff outer")) %>% 
        mutate(fdr = ifelse(is.na(fdr), 0, fdr), 
               frac_correct_rejections = ifelse(is.na(frac_correct_rejections), 0, frac_correct_rejections), 
               n_outcomes_implicated = ifelse(is.na(n_outcomes_implicated), 0, n_outcomes_implicated)) %>% # FDP is missing if power = 0, so also put to 0
        summarise(power = mean(power),  
                  fdr = mean(fdr), 
                  gamma = mean(gamma), 
                  n_outcomes_implicated = mean(n_outcomes_implicated), 
                  frac_correct_rejections = mean(frac_correct_rejections)) %>% 
        mutate(sparsity_plot = paste0("sparsity = ", s))
      
      power_fdp <- rbind(power_fdp, snr_res_averaged)
      
    }
    
    # check gamma 
    as.data.frame(power_fdp) %>% 
      ggplot(aes(x = np_ratio, y = gamma, col = method)) + #
      geom_line(size = 1.6) +
      labs(x = "n/p", y = "gamma", color = "Method") + 
      theme_minimal() + 
      theme(text = element_text(size=25))
    
    
    across_resolution = power_fdp %>% 
      dplyr::select(method, np_ratio, power, fdr, frac_correct_rejections, n_outcomes_implicated, sparsity_plot)
    
    ########## KNOCKOFF BY RESOLUTION
    power_fdp <- data.frame()
    
    my_files <- list.files(mydir)
    my_files <- my_files[str_detect(my_files, "hp_union_by_res")]
    my_files <- my_files[str_detect(my_files, paste0("nrep_", nrep, "_"))]
    my_files <- my_files[str_detect(my_files, paste0("sp_", s))]
    my_files <- my_files[str_detect(my_files, paste0("ol_", ol))]
    my_files <- my_files[str_detect(my_files, paste0("ntg", ntg))]
    my_files <- my_files[str_detect(my_files, paste0("dpo", po))]
    my_files <- my_files[str_detect(my_files, paste0("p", p))]
    
    ldf <- lapply(my_files, read.csv)
    
    for(l in 1:length(ldf)) {
      
      
      snr_res_averaged <- ldf[[l]]  %>% 
        dplyr::filter(np_ratio <= np_ratio_max) %>%
        group_by(method, resolution, np_ratio) %>% 
        mutate(fdr = ifelse(is.na(fdp), 0, fdp), 
               frac_correct_rejections = ifelse(is.na(frac_correct), 0, frac_correct), 
               n_outcomes_implicated = ifelse(is.na(n_outcomes_implicated), 0, n_outcomes_implicated)) %>% # FDP is missing if power = 0, so also put to 0
        summarise(power = mean(power),  
                  fdr = mean(fdr), 
                  gamma = mean(gamma), 
                  n_outcomes_implicated = mean(n_outcomes_implicated), 
                  frac_correct_rejections = mean(frac_correct)) %>% 
        mutate(sparsity_plot = paste0("sparsity = ", s))
      
      
      power_fdp <- rbind(power_fdp, snr_res_averaged)
      
    }
    
    
    knockoff_by_res <- power_fdp %>% 
      mutate(method = paste0(method, " res.",resolution)) %>% 
      dplyr::filter(method != "ebh knockoff res.3") %>%
      dplyr::select(method, np_ratio,power, fdr, frac_correct_rejections, n_outcomes_implicated, sparsity_plot)
    knockoff_by_res <- knockoff_by_res[, c(2, 3, 4, 5, 6, 7, 8)]
    
    
    
    combined <- rbind(across_resolution, knockoff_by_res)
    
  
    
    all_res <- rbind(all_res, combined)
    
    
    
    
    
    
  }
  

  color_values <- c("steelblue1", "mediumpurple3", "mediumpurple3", 
                    "orange4", "orange4", "darkorange3",  "darkorange3","darkgoldenrod1")
  
  
  linetype = c("solid", "dotdash", "solid", "dotdash","solid", "dotdash", "solid",  "solid")
  
  all_res <- all_res %>% 
    mutate(method = factor(method, levels = c("e-FBH", 
                                              "knockoff outer", 
                                              "ebh knockoff outer",
                                              "knockoff res.1", 
                                              "ebh knockoff res.1",
                                              "knockoff res.2", 
                                              "ebh knockoff res.2",
                                              "knockoff res.3"))) %>%
    mutate(method = recode(method, `e-FBH` = "KeLP",
                           `knockoff outer` = "Knockoff outer", 
                           `ebh knockoff outer` = "e-BH knockoff outer", 
                           `ebh knockoff res.1` = "e-BH knockoff leaves", 
                           `ebh knockoff res.2` = "e-BH knockoff internal",
                           `knockoff res.1` = "Knockoff leaves",
                           `knockoff res.2` = "Knockoff internal node",
                           `knockoff res.3` = "Knockoff root node")) 
  
  as.data.frame(all_res) %>% 
    ggplot(aes(x = np_ratio, y = frac_correct_rejections, col = method, linetype = method)) + facet_grid(~sparsity_plot) +
    geom_line(size = 1.6) +
    labs(x = "n/p", y = "Power", color = "Method", linetype = "Method") + 
    theme_light() + 
    scale_color_manual(values=color_values, guide = guide_legend(nrow = 3,  byrow = TRUE)) +
    scale_linetype_manual(values = linetype) + 
    theme(text = element_text(size=25)) + 
    theme(strip.background = element_blank(), strip.text.x = element_blank()) -> np_ratio_frac_correct
  
  as.data.frame(all_res) %>% 
    ggplot(aes(x = np_ratio, y = fdr, col = method, linetype = method)) + facet_grid(~sparsity_plot) +
    geom_line(size = 1.6) +
    geom_hline(yintercept = alpha, linetype = "dashed", col = "plum") + 
    labs(x = "n/p", y = "FDR", color = "Method", linetype = "Method") + 
    theme_light() + 
    scale_color_manual(values=color_values, guide = guide_legend(nrow = 3,  byrow = TRUE)) +
    scale_linetype_manual(values = linetype) + 
    theme(text = element_text(size=25)) + 
    theme(strip.text = element_text(colour = 'black', size = 25)) -> np_ratio_fdr

  as.data.frame(all_res) %>% 
    ggplot(aes(x = frac_correct_rejections, y = sqrt(n_outcomes_implicated), col = method, linetype = method)) + facet_grid(~sparsity_plot) +
    geom_line(size = 1.6) +
    labs(x = "Power", y = "Sqrt(Outcomes implicated)", color = "Method", linetype = "Method") + 
    theme_light() + 
    scale_color_manual(values=color_values, guide = guide_legend(nrow = 3,  byrow = TRUE)) +
    scale_linetype_manual(values = linetype) + 
    theme(text = element_text(size=25)) + 
    theme(strip.background = element_blank(), strip.text.x = element_blank()) -> frac_correct_outcomes_implicated
  
 
  ggarrange(np_ratio_fdr, np_ratio_frac_correct, frac_correct_outcomes_implicated, common.legend = TRUE, legend="bottom", nrow = 3) -> combined_plot
  combined_plot
  ggsave(paste0(my_out_dir,"figures/hierarchical_combined_plot_", nrep, "_a_", alpha, "_ol_", ol, "_p_",p,"_po_",po,"_ntg_",ntg,".png"), combined_plot, 
         width = 15, height = 12)
  
  
  
}    

}
  
}
