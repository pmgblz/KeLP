


# set paths and simulation parameters

if(dir.exists("~")){
  mydir <- "~"
  setwd(mydir)
  
  my_out_dir <- "~" 
}

settings <- c("0.025", "0.05","0.1") 

gammas <- "tuned"

amp_max = 30

chromosome <- 21


nrep <- 25 
alpha <- 0.2  

tg = 2

all_res <- data.frame()
all_snrs <- data.frame()

for(s in settings) {

  
  # KNOCKOFF
  my_files <- list.files(mydir)
  my_files <- my_files[str_detect(my_files, "UKB_kelp_res")]
  my_files <- my_files[str_detect(my_files, paste0("_c", chromosome))]
  my_files <- my_files[str_detect(my_files, "res_by")]
  my_files <-  my_files[str_detect(my_files, paste0("nrep_", nrep))]
  my_files <-  my_files[str_detect(my_files, paste0("sp", s))]
  my_files <-  my_files[str_detect(my_files, paste0("tg", tg))]
  
  
  ldf <- lapply(my_files, read.csv)
  
  power_fdp <- data.frame()
  
  
  for(l in 1:length(ldf)) {
    
    
    snr_res_averaged <- ldf[[l]]  %>% 
      #filter(Resolution == 1) %>% # consider the first layer only
      group_by(method, Resolution, amp, Q) %>% #
      dplyr::filter(amp  <= amp_max) %>%
      mutate(fdp = ifelse(is.na(fdp), 0, fdp), 
             power = ifelse(is.na(power), 0, power), 
             snps_implicated = ifelse(is.na(snps_implicated), 0, snps_implicated), 
             frac_correct_rejections = ifelse(is.na(frac_correct_rejections), 0, frac_correct_rejections)) %>% # FDP is missing if power = 0, so also put to 0
      summarise(power = mean(power), 
                fdp = mean(fdp), 
                snps_implicated = mean(snps_implicated), 
                frac_correct_rejections = mean(frac_correct_rejections))
    
    power_fdp <- rbind(power_fdp, snr_res_averaged)
    
  }
  
  
knockoff_res <- power_fdp %>% dplyr::mutate(method = paste0("Knockoff res. ", Resolution)) 
  knockoff_res <- knockoff_res[, c(1, 3, 6, 7, 8)]
  # KELP
  
  
  my_files <- list.files(mydir)
  my_files <- my_files[str_detect(my_files, "UKB_kelp_res")]
  my_files <- my_files[str_detect(my_files, "res_across")]
  my_files <-  my_files[str_detect(my_files, paste0("nrep_", nrep))]
  my_files <-  my_files[str_detect(my_files, paste0("sp", s))]
  my_files <-  my_files[str_detect(my_files, paste0("tg", tg))]
  
  ldf <- lapply(my_files, read.csv)
  
  power_fdp <- data.frame()
  
  
  for(l in 1:length(ldf)) {
    
    
    snr_res_averaged <- ldf[[l]]  %>% 
      #filter(Resolution == 1) %>% # consider the first layer only
      dplyr::filter(amp  <= amp_max) %>%
      group_by(method, amp, Q) %>% #
      mutate(fdp = ifelse(is.na(fdp), 0, fdp), 
             power = ifelse(is.na(power), 0, power), 
             snps_implicated = ifelse(is.na(snps_implicated), 0, snps_implicated), 
             frac_correct_rejections = ifelse(is.na(frac_correct_rejections), 0, frac_correct_rejections)) %>% # FDP is missing if power = 0, so also put to 0
      summarise(power = mean(power), 
                fdp = mean(fdp), 
                snps_implicated = mean(snps_implicated), 
                frac_correct_rejections = mean(frac_correct_rejections))
    
    power_fdp <- rbind(power_fdp, snr_res_averaged)
    
  }
  
  
  
  kelp_res <- power_fdp[, c(1, 2, 5, 6, 7)]
  
  
  
  knockoff_kelp_combined <- rbind(knockoff_res, kelp_res)  
  
  knockoff_kelp_combined <- knockoff_kelp_combined %>% 
    mutate(method = factor(method, levels = c("kelp", 
                                              "knockoffouter", 
                                              "ebhknockout",
                                              "Knockoff res. 0",
                                              "Knockoff res. 2",
                                              "Knockoff res. 4",
                                              "Knockoff res. 6"))) %>%
    mutate(method = recode(method, kelp = "KeLP", knockoffouter = "Knockoff outer", ebhknockout = "e-BH knockoff outer", 
                           `Knockoff res. 0` = "Knockoff single-SNP", `Knockoff res. 2` = "Knockoff 20kb", 
                           `Knockoff res. 4` = "Knockoff 81kb", `Knockoff res. 6` = "Knockoff 425kb")) %>% 
    mutate(sparsity_plot = paste0("sparsity = ", s))
  
  
  all_res <- rbind(all_res, knockoff_kelp_combined)
  
  
  # get SNR
  # KNOCKOFF
  my_files <- list.files(mydir)
  my_files <- my_files[str_detect(my_files, "UKB_")]
  my_files <- my_files[str_detect(my_files, paste0("_c", chromosome))]
  my_files <- my_files[str_detect(my_files, "_snr_")]
  my_files <-  my_files[str_detect(my_files, paste0("nrep_", nrep))]
  my_files <-  my_files[str_detect(my_files, paste0("sp", s))]
  my_files <-  my_files[str_detect(my_files, paste0("tg", tg))]
  
  
  ldf <- lapply(my_files, read.csv)
  
  
  
  
  for(l in 1:length(ldf)) {
    
    
    snr_res_averaged <- ldf[[l]]  %>% 
      #filter(Resolution == 1) %>% # consider the first layer only
      group_by(amp) %>% #
      summarise(mean_snr = mean(my_snr)) %>% 
      mutate(sparsity_plot = paste0("sparsity = ", s)) %>% 
      mutate(mean_snr_alt = mean_snr / (1 + mean_snr))
    
    all_snrs <- rbind(all_snrs, snr_res_averaged)
    
  }
  
  
  
}


color_values <- c("steelblue1","mediumpurple3", "mediumpurple3","orange4", "darkorange3","darkgoldenrod1", "gold1")

linetype = c("solid", "dotdash", "solid", "solid", "solid", "solid", "solid")

  as.data.frame(all_res) %>% 
    ggplot(aes(x = frac_correct_rejections, y = sqrt(snps_implicated), color = method, linetype = method)) + 
    facet_grid(~sparsity_plot) +
    geom_line(size = 1.8) + 
    labs(x = "Power", y = "sqrt(SNPs implicated)", col = "Method", linetype = "Method") + 
    theme_light()  +
    theme(text = element_text(size=30))  +
    scale_color_manual(values=color_values, guide = guide_legend(nrow = 3,  byrow = TRUE)) + 
    scale_linetype_manual(values = linetype) +
    theme(strip.background = element_blank(), strip.text.x = element_blank())  -> frac_snps

  as.data.frame(all_res) %>% 
    ggplot(aes(x = amp, y = sqrt(snps_implicated), color = method, linetype = method)) + 
    facet_grid(~sparsity_plot) +
    geom_line(size = 1.8) +
    labs(x = "Signal Amplitude", y = "sqrt(SNPs implicated)", col = "Method", linetype = "Method") + 
    theme_light()  +
    theme(text = element_text(size=30))  +
    scale_color_manual(values=color_values, guide = guide_legend(nrow = 3,  byrow = TRUE))  +
    scale_linetype_manual(values = linetype) +
    theme(strip.background = element_blank(), strip.text.x = element_blank())
  
  as.data.frame(all_res) %>% 
    ggplot(aes(x = amp, y = fdp, color = method, linetype = method))  +
    facet_grid(~sparsity_plot) +
    coord_cartesian(ylim = c(0, max(all_res$fdp))) +
    geom_hline(yintercept = alpha, linetype = "dashed", col = "black") +
    geom_line(size = 1.8) +
    labs(x = "Signal amplitude", y = "FDR", col = "Method", linetype = "Method") + 
    theme_light()  +
    theme(text = element_text(size=30))   +
    scale_color_manual(values=color_values, guide = guide_legend(nrow = 3,  byrow = TRUE)) + 
    scale_linetype_manual(values = linetype) +
    theme(strip.text = element_text(colour = 'black', size = 30))  -> amp_fdp
  
  as.data.frame(all_res) %>% 
    ggplot(aes(x = amp, y = frac_correct_rejections, color = method, linetype = method)) + 
    facet_grid(~sparsity_plot) +
    geom_line(size = 1.8) +
    labs(x = "Signal amplitude", y = "Power", col = "Method", linetype = "Method") + 
    theme_light()  +
    theme(text = element_text(size=30))  +
    scale_color_manual(values=color_values, guide = guide_legend(nrow = 3,  byrow = TRUE)) + 
    scale_linetype_manual(values = linetype) +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) -> amp_frac_correct
  
  
  as.data.frame(all_snrs) %>% 
    ggplot(aes(x = amp, y = mean_snr / (mean_snr + 1))) + 
    facet_grid(~sparsity_plot) +
    geom_line(size = 1.8) +
    labs(x = "Signal amplitude", y = "Var(f(x))/Var(Y)") + 
    theme_light()  +
    theme(text = element_text(size=30))  +
    coord_cartesian(ylim = c(0, 1)) +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) -> snr_plot
  
  
  # power_by_sparsity
  ggarrange(amp_fdp, amp_frac_correct, frac_snps, snr_plot, common.legend = TRUE, legend="bottom", nrow = 4) -> combined_plot
  combined_plot
  ggsave(paste0(my_out_dir,"figures/UKB_combined_plot_nrep",nrep,"_chr_", chromosome, "_tg", tg ,".png"), combined_plot, 
         width = 18, height = 16)
  



