
###################################################################### 

# Summarize kelp simulation results and produce power / fdp plots #

###################################################################### 

rm(list = ls())


if(dir.exists("")){
  mydir <- ""
  setwd(mydir)
  
  my_out_dir <- "" 
}


library(ggpubr)


nrep <- 100
alpha <- 0.2 


bsd <- c("0.2") 

sparsity <- c(0.025, 0.05, 0.1) 


rhos <- c(0.8)

p <- 1000

np_max <- 3

for(r in rhos) {

for(gtt in c("1.0")) {
  
    for(b in bsd) {
      
      for(c in c(2)) {
        
        
        
        all_frac_rejected_combined = c()
        
        for(s in sparsity) {
          
    
          
          
          print(b)
          print(c)
          print(s)
          
          
          ########### GET AVERAGE SNR ##############
          
          # GET SNR 
          my_files <- list.files(mydir)
          my_files <- my_files[str_detect(my_files, "np")]
          my_files <- my_files[str_detect(my_files, paste0("p", p))]  
          my_files <- my_files[str_detect(my_files, paste0("r", r))]  
          my_files <- my_files[str_detect(my_files, "_snr_")]
          my_files <- my_files[str_detect(my_files, paste0("_ncuts_", c))]
          my_files <- my_files[str_detect(my_files, paste0("_gtt", gtt))]
          my_files <- my_files[str_detect(my_files, paste0("s", s))]
          my_files <- my_files[str_detect(my_files, paste0("sdb", b))]   
          
          my_files <-  my_files[str_detect(my_files, paste0("nrep_", nrep))]

          
          my_files <-  my_files[str_detect(my_files, "Q_1")]
          
          
          
          
          Sys.sleep(1)  
          ldf <- lapply(my_files, read.csv)
          Sys.sleep(1)  
          all_snr <- data.frame()
          
          
          for(l in 1:length(ldf)) {
            
            
            snr_res_averaged <- mean(ldf[[l]]$my_snr)
            
            all_snr <- rbind(all_snr, snr_res_averaged)
            
          }
          
          avg_snr <- mean(all_snr[, 1])
          print("SNR IS: ")
          print(avg_snr)
          
          
          
          ########### COMPARE FRAC REJECTED AND SNPS IMPLICATED #############
          ############ ONLY FOR Q = 1 ###############
          
          # FOR THE KNOCKOFF, GET SNPs implicated and fraction correct for each level of resolution
          
          Sys.sleep(1)  
          
          my_files <- list.files(mydir)
          my_files <- my_files[str_detect(my_files, "np")]
          my_files <- my_files[str_detect(my_files, paste0("p", p))]  
          my_files <- my_files[str_detect(my_files, paste0("r", r))]  
          my_files <- my_files[!str_detect(my_files, "res_across")]
          my_files <- my_files[!str_detect(my_files, "_snr_")]
          my_files <- my_files[str_detect(my_files, paste0("_ncuts_", c))]
          my_files <- my_files[str_detect(my_files, paste0("_gtt", gtt))]
          my_files <- my_files[str_detect(my_files, paste0("s", s))]
          my_files <- my_files[str_detect(my_files, paste0("sdb", b))]   

          
          my_files <-  my_files[str_detect(my_files, paste0("nrep_", nrep))]

          my_files <-  my_files[str_detect(my_files, "Q_1")]

          Sys.sleep(1)  
          ldf <- lapply(my_files, read.csv)
          Sys.sleep(1)  
          power_fdp <- data.frame()
          
          
          for(l in 1:length(ldf)) {
            
            
            snr_res_averaged <- ldf[[l]]  %>% 
              group_by(method, Resolution, np_ratio, Q) %>% #
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
          
          
          knockoff_snps_implicated_frac_correct_by_res <- power_fdp %>% dplyr::filter(method == "knockoff") %>% 
            mutate(method = paste0("Knockoff res. ", Resolution)) %>% 
            dplyr::select(method, frac_correct_rejections, snps_implicated, np_ratio, Resolution, fdp, power)
          knockoff_snps_implicated_frac_correct_by_res <- knockoff_snps_implicated_frac_correct_by_res[, c(1, 2, 3, 4, 6, 7)]
          colnames(knockoff_snps_implicated_frac_correct_by_res) <- c("method", "frac_correct_rejections", "snps_implicated", "np_ratio", "fdp", "power") 
   
          # now for kelp
          my_files <- list.files(mydir)
          my_files <- my_files[str_detect(my_files, "np")]
          my_files <- my_files[str_detect(my_files, paste0("p", p))]  
          my_files <- my_files[str_detect(my_files, paste0("r", r))]  
          my_files <- my_files[str_detect(my_files, "res_across")]
          my_files <- my_files[!str_detect(my_files, "_snr_")]
          my_files <- my_files[str_detect(my_files, paste0("_ncuts_", c))]
          my_files <- my_files[str_detect(my_files, paste0("_gtt", gtt))]
          my_files <- my_files[str_detect(my_files, paste0("s", s))]
          my_files <- my_files[str_detect(my_files, paste0("sdb", b))]   
          my_files <-  my_files[str_detect(my_files, paste0("nrep_", nrep))]
          
          ldf <- lapply(my_files, read.csv)
          
          Sys.sleep(2)  
          
          power_fdp <- data.frame()
          
          
          for(l in 1:length(ldf)) {
            
            
            snr_res_averaged <- ldf[[l]]  %>%
              group_by(method, np_ratio, Q) %>% 
              mutate(fdp = ifelse(is.na(fdp), 0, fdp), 
                     power = ifelse(is.na(power), 0, power)) %>% # FDP is missing if power = 0, so also put to 0
              summarise(power = mean(power), 
                        fdp = mean(fdp), 
                        snps_implicated = mean(snps_implicated, na.rm = TRUE), 
                        frac_correct_rejections = mean(frac_correct_rejections, na.rm = TRUE))
            
            power_fdp <- rbind(power_fdp, snr_res_averaged)
            
          }
          
          kelp_snps_implicated_frac_correct_by_res <- power_fdp %>% 
            dplyr::select(method, frac_correct_rejections, snps_implicated, np_ratio, fdp, power)
          
          kelp_snps_implicated_frac_correct_by_res <- kelp_snps_implicated_frac_correct_by_res[, c(1, 2, 3, 4, 5, 6)]
          
          
          frac_rejected_combined <- rbind(kelp_snps_implicated_frac_correct_by_res, knockoff_snps_implicated_frac_correct_by_res) 
          frac_rejected_combined <- frac_rejected_combined %>% mutate(sparsity = s) %>% mutate(avg_snr_rounded = round(avg_snr, 1))
      
          
          
          all_frac_rejected_combined <- rbind(all_frac_rejected_combined, frac_rejected_combined)
          
          
          
        }
        
        
        
 
       
        
        
        
        if(c == 3) {
          linetype = c("solid", "dotdash", "solid", "solid", "solid", "solid")
          color_values <- c("cornflowerblue","mediumpurple3", "mediumpurple3","darkgoldenrod", "darkgoldenrod1","burlywood1")
          all_frac_rejected_combined <- all_frac_rejected_combined %>% 
            mutate(method = factor(method, levels = c("kelp", 
                                                      "knockoff_outer", 
                                                      "knockoff_outer_ebh",
                                                      "Knockoff res. 1",
                                                      "Knockoff res. 2",
                                                      "Knockoff res. 3"))) %>%
            mutate(method = recode(method, kelp = "KeLP", knockoff_outer = "Knockoff outer", knockoff_outer_ebh = "e-BH knockoff outer")) 
        } else if(c == 4) {
          linetype = c("solid", "dotdash", "solid", "solid", "solid", "solid", "solid")
          color_values <- c("cornflowerblue","mediumpurple3", "mediumpurple3","darkgoldenrod", "darkgoldenrod1","burlywood1", "cornsilk")
          all_frac_rejected_combined <- all_frac_rejected_combined %>% 
            mutate(method = factor(method, levels = c("kelp", 
                                                      "knockoff_outer", 
                                                      "knockoff_outer_ebh",
                                                      "Knockoff res. 1",
                                                      "Knockoff res. 2",
                                                      "Knockoff res. 3", 
                                                      "Knockoff res. 4"))) %>%
            mutate(method = recode(method, kelp = "KeLP", knockoff_outer = "Knockoff outer", knockoff_outer_ebh = "e-BH knockoff outer")) 
        } else {
          #color_values <- c("cornflowerblue","mediumpurple3", "mediumpurple3","darkgoldenrod", "darkgoldenrod1")
          color_values <- c("steelblue1","mediumpurple3", "mediumpurple3","orange4", "gold1")
          linetype = c("solid", "dotdash", "solid", "solid", "solid")
          all_frac_rejected_combined <- all_frac_rejected_combined %>% 
            mutate(method = factor(method, levels = c("kelp", 
                                                      "knockoff_outer", 
                                                      "knockoff_outer_ebh",
                                                      "Knockoff res. 1",
                                                      "Knockoff res. 2"))) %>%
            mutate(method = recode(method, kelp = "KeLP", 
                                   knockoff_outer = "Knockoff outer",
                                   knockoff_outer_ebh = "e-BH knockoff outer", 
                                   `Knockoff res. 1` = "Knockoff individual",
                                   `Knockoff res. 2` = "Knockoff group")) 
        }
        
        
        all_frac_rejected_combined %>% 
          mutate(sparsity_plot = paste0("sparsity = ", sparsity)) -> all_frac_rejected_combined
        
        
        all_frac_rejected_combined %>% 
          dplyr::filter(np_ratio <= np_max) %>%
          ggplot(aes(x = np_ratio, y = frac_correct_rejections, col = method, linetype = method)) + facet_grid(~sparsity_plot) + 
          geom_line(size = 1.5) + 
          theme_light()  +
          coord_cartesian(ylim = c(0, 1)) +
          theme(panel.spacing = unit(2, "lines")) +
          labs(x = "n/p", y = "Power", color = "Method", linetype = "Method") + 
          scale_color_manual(values=color_values, guide = guide_legend(nrow = 2,  byrow = TRUE)) + 
          scale_linetype_manual(values = linetype) +
          theme(text = element_text(size=25))  +
          theme(strip.background = element_blank(), strip.text.x = element_blank())  -> frac_correct_by_sparsity
        
        
        all_frac_rejected_combined %>% 
          dplyr::filter(np_ratio <= np_max) %>%
          ggplot(aes(x = frac_correct_rejections, y = snps_implicated, col = method, linetype = method)) + facet_grid(~sparsity_plot) + 
          geom_line(size = 1.5) + 
          theme_light()  +
          theme(text = element_text(size=25)) +
          theme(panel.spacing = unit(2, "lines")) +
          labs(x = "Power", y = "sqrt(SNPs implicated)", color = "Method", linetype = "Method") + 
          scale_color_manual(values=color_values, guide = guide_legend(nrow = 2,  byrow = TRUE))  +
          scale_linetype_manual(values = linetype) +
          theme(text = element_text(size=25))  +
          theme(strip.background = element_blank(), strip.text.x = element_blank()) -> snps_correct_by_sparsity
        
        all_frac_rejected_combined %>% 
          dplyr::filter(np_ratio <= np_max) %>%
          ggplot(aes(x = np_ratio, y = fdp, col = method, linetype = method)) + facet_grid(~sparsity_plot) + 
          geom_line(size = 1.5) + 
          theme_light()  +
          theme(text = element_text(size=25)) +
          theme(panel.spacing = unit(2, "lines")) +
          geom_hline(yintercept = 0.2, linetype = "dashed", col = "black") + 
          labs(x = "n/p", y = "FDR", color = "Method", linetype = "Method") + 
          scale_color_manual(values=color_values, guide = guide_legend(nrow = 2,  byrow = TRUE)) + 
          scale_linetype_manual(values = linetype) +
          theme(text = element_text(size=25))  +
          theme(strip.text = element_text(colour = 'black', size = 25)) -> fdr_by_sparsity

        # power_by_sparsity
        ggarrange(fdr_by_sparsity, frac_correct_by_sparsity, snps_correct_by_sparsity, common.legend = TRUE, legend="bottom", nrow = 3) -> combined_plot
        combined_plot
        ggsave(paste0(my_out_dir,"figures/block_np_combined_plot_", nrep, "_a_", alpha,"_nc_",c, "_bsd_",b , "_p_", p, "_gtt_",gtt ,"_rho_",r,".png"), combined_plot, 
               width = 13, height = 10)
        
        
        print("PLOT PRODUCED DONE")
        
        
        

        
        
      }  
      
    }  
    
}

}
