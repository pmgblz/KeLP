

###################################################################### 

# Summarize e-MLKF simulation results and produce power / fdp plots #

###################################################################### 



if(dir.exists("")){
  mydir <- ""
  setwd(mydir)
  
  my_out_dir <- "" 
}


nrep <- 100
alpha <- 0.2    
amp_max = 12
k <- c(10, 20)
M <- 2
ck_mult_ending <- c("mult")
gamma_ending <- c("tuned") 


for(nz in k) {
  
  
  for(ck in ck_mult_ending) {
    
    if(ck == "not_mult") {
      filename_ck_ending <- "notmult"
    } else {
      filename_ck_ending <- "mult"
    }
    
    
    for(ge in gamma_ending) {
      
   
    
    my_files <- list.files(mydir)
    my_files <- my_files[str_detect(my_files, paste0("_", M, "M"))]
    my_files <-  my_files[str_detect(my_files, paste0("nrep_", nrep))]
    my_files <-  my_files[str_detect(my_files, paste0("nnbg_" ,nz))]
    my_files <-  my_files[str_detect(my_files, paste0("_ck_" ,ck))]
    my_files <-  my_files[str_detect(my_files, paste0("gamma_" ,ge))]
    
    
    ldf <- lapply(my_files, read.csv)
    
    power_fdp <- data.frame()
    
    
    for(l in 1:length(ldf)) {
      
      
      snr_res_averaged <- ldf[[l]]  %>% 
        group_by(Method, Layer, Amplitude) %>% #
        filter(Amplitude <= amp_max) %>%
        mutate(fdp = ifelse(is.na(FDP), 0, FDP))  %>% # FDP is missing if power = 0, so also put to 0
        summarise(power = mean(Power), 
                  fdp = mean(fdp), 
                  gamma = mean(Gamma, na.rm = TRUE)) %>% 
        mutate(Method = ifelse(Method == "EMLKF", "e-MKF", 
                               ifelse(Method == "MLKF", "MKF", Method)))
      
      power_fdp <- rbind(power_fdp, snr_res_averaged)
      
    }
    
    
    
    as.data.frame(power_fdp) %>% 
      ggplot(aes(x = Amplitude, y = power, col = Method, linetype = as.factor(Layer))) + #
      geom_line(size = 1.6) +
      labs(x = "Signal amplitude", y = "Power", color = "Method", linetype = "Layer" ) + 
      theme_minimal() + 
      theme(text = element_text(size=25))  + 
      scale_color_manual(values=c("steelblue1", "orange4"))
    
    ggsave(paste0(my_out_dir,"figures/emlkf_power_nrep", nrep, "_a_", alpha, "ck", filename_ck_ending, "_k_",nz ,"_M", M, "_gamma_", ge,".png"), width = 7, height = 6)
    
    
    alpha <- 0.2
    
    
    ylim_fdr <- 0 
    if(max(power_fdp$fdp) <= alpha) {
      ylim_fdr <- alpha + 0.05
    } else {
      ylim_fdr <- max(power_fdp$fdp) + 0.05
    }
    
    
    as.data.frame(power_fdp) %>% 
      ggplot(aes(x = Amplitude, y = fdp, color = Method , linetype = as.factor(Layer))) + # , 
      geom_line(size = 1.6)+ 
      labs(x = "Signal amplitude", y = "FDR", color = "Method" , linetype = "Layer") + #
      theme_minimal()  +  coord_cartesian(ylim = c(0, ylim_fdr)) +
      geom_hline(yintercept = alpha, linetype = "dashed", col = "plum") + 
      theme(text = element_text(size=25)) + 
      scale_color_manual(values=c("steelblue1", "orange4"))
    
    ggsave(paste0(my_out_dir,"figures/emlkf_fdr_nrep", nrep, "_a_", alpha, "ck", filename_ck_ending, "_k_",nz ,"_M", M, "_gamma_",ge,".png"), width = 7, height = 6)
    
    rm(ldf)
    rm(power_fdp)
    rm(snr_res_averaged)
    rm(my_files)
    
    }
    
  }
  
 
}





