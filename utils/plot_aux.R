
                      # PLOT FUNCTIONS FOR UKB #

# This file contains function and parameters for plotting UKB rejections 
# code based on M. Sesia, S. Bates, E. Candès, J. Marchini, C. Sabatti Proceedings of the National Academy of Sciences, 2021; doi:10.1073/pnas.2105841118 
# https://github.com/msesia/knockoffgwas/blob/master/visualization/utils_plotting.R


# some baseline blotting parameters from https://github.com/msesia/knockoffgwas/blob/master/visualization/utils_plotting.R
font.size <- 18
title.font.size <- 18
axis.font.size <- 18
legend.font.size <- 10
dot.size <- 2.5

bp.labeler_simulation <- function(x) {
  round(x,0)
}

bp.labeler_UKB <- function(x) {
  round(x*1e-6,3)
}

# Function to plot a chicago plot for simulation results
  # INPUT: 
  # Discoveries: data frame containing information on the discoveries made by a method 
  # multiple_u: indicator for whether multiple u where used in partial conjunction
  # OUTPUT: 
  # plot_chicago: chicago plot

plot_chicago_simulation <- function(Discoveries, multiple_u = FALSE) {
  
  window.left = min(Discoveries$BP.min)
  window.right = max(Discoveries$BP.max)
  
  # Extract knockoff discoveries within this window
  if(!is.null(Discoveries)) {
    Knockoffs.window <- Discoveries %>%
      filter(BP.min<=window.right, BP.max>=window.left)
  } else {
    Knockoffs.window <- tibble()
  }
  
  # Plot knockoff discoveries
  resolution.list <- c("res1", "res2", "res3") %>% rev
  resolution.heights <- seq(length(resolution.list))
  names(resolution.heights) <- resolution.list
  resolution.labels <- resolution.list
  
  if(nrow(Knockoffs.window)>0) {
    
    
    if(multiple_u) {
      plot_chicago <- Knockoffs.window %>%
        mutate(Resolution=as.character(resolution)) %>%
        mutate(Resolution = paste0("res", Resolution)) %>%
        mutate(Height=resolution.heights[Resolution]) %>%
        mutate(`Part. Conj. u` = as.factor(u)) %>%
        ggplot() +
        geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5, fill = `Part. Conj. u`, 
                      color = `Part. Conj. u`)) + 
        scale_fill_manual(values=c("plum4","#E69F00", "#56B4E9")) + 
        scale_color_manual(values=c("plum4","#E69F00", "#56B4E9"))
    } else {
      plot_chicago <- Knockoffs.window %>%
        mutate(Resolution=as.character(resolution)) %>%
        mutate(Resolution = paste0("res", Resolution)) %>%
        mutate(Height=resolution.heights[Resolution]) %>%
        ggplot() +
        geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5), color = "black")
    }
   
  } else {
    plot_chicago <- ggplot(tibble()) + geom_blank()
  }
  
  plot_chicago <- plot_chicago +
    labs(x = "", y = "Resolution") +
    coord_cartesian(xlim = c(window.left,window.right)) +
    scale_x_continuous(expand=c(0.01,0.01), labels=bp.labeler_simulation) +
    scale_y_continuous(limits=c(0.5,max(resolution.heights)+0.5),
                       labels=resolution.labels, breaks=resolution.heights) +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          panel.border=element_blank(),
          panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
          panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray"),
          text = element_text(size=font.size),
          axis.title.y = element_text(size=title.font.size),
          plot.title = element_text(size=title.font.size),
          legend.text = element_text(size=legend.font.size),
          legend.title = element_text(size=legend.font.size),
          legend.key.height = unit(0.75,"line")
    )
  
  return(plot_chicago)
}

########### PLOT CHICAGO FOR UKB DATA ####################### 

# two functions: general chicago plot, chicago plot including implicitly rejected groups 

# Function to plot a chicago plot for UKB results
# INPUT: 
  # window.chr: which chromosome to plot 
  # window.left: BP starting position for plot
  # window.right: BP end position for plot
  # Discoveries: data frame containing information on the discoveries made by a method 
  # evals: indicator for whether to fill based on evalues or evalues without multiplier (fractions)
# OUTPUT: 
  # plot_chicago: chicago plot


plot_chicago <- function(window.chr, window.left, window.right, Discoveries, evals = FALSE, fill_in = TRUE) {
  
  # Extract knockoff discoveries within this window
  if(!is.null(Discoveries)) {
    Knockoffs.window <- Discoveries %>%
      dplyr::filter(CHR==window.chr, BP.min<=window.right, BP.max>=window.left)
    #cat(sprintf("There are %d knockoff discoveries within this window.\n", nrow(Knockoffs.window)))
  } else {
    Knockoffs.window <- tibble()
  }
  
  # Plot knockoff discoveries
  resolution.list <- c("res0", "res1", "res2", "res3", "res4", "res5", "res6") %>% rev
  resolution.heights <- seq(length(resolution.list))
  names(resolution.heights) <- resolution.list
  #resolution.labels <- paste(parse_number(resolution.list), "\\%", sep="")
 # resolution.labels <- resolution.list
  resolution.labels <- c("single-SNP", "3 kb", "20 kb", "41 kb", "81 kb", "208 kb", "425 kb") %>% rev
  
  #resolution.labels <- c("single-SNP", "3", "20", "41", "81", "208", "425") %>% rev
  
  if(nrow(Knockoffs.window)>0) {
    
    if(evals & fill_in) {
      
      plot_chicago <- Knockoffs.window %>%
        mutate(Resolution=as.character(Resolution)) %>%
        mutate(Resolution = paste0("res", Resolution)) %>%
        mutate(Height=resolution.heights[Resolution]) %>%
        ggplot() +
        geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5, fill=evals - mean(evals)), #, 
                  color="black")  +
        scale_fill_gradient(name="Evalue",
                            #low="red1", high="dodgerblue1",
                            low="gray98", high="gray30",
                            limits=c(0,  round(quantile(Knockoffs.window$evals, 0.9)[[1]])),
                            labels = c("0", "q50", "q75", "q90"),
                            breaks=c(0,round(quantile(Knockoffs.window$evals, 0.5)[[1]]),
                                     round(quantile(Knockoffs.window$evals, 0.75)[[1]]), 
                                     round(quantile(Knockoffs.window$evals, 0.9)[[1]])),
                            space="Lab", na.value="gray", guide="colourbar") # ,labels=function(x){x}
      
    }
    
    
    if(!evals & fill_in) {
      
      plot_chicago <- Knockoffs.window %>%
        mutate(Resolution=as.character(Resolution)) %>%
        mutate(Resolution = paste0("res", Resolution)) %>%
        mutate(Height=resolution.heights[Resolution]) %>%
        ggplot() +
        geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5, fill=fracs - mean(fracs)), #, 
                  color="black")  +
        scale_fill_gradient(name="E-value",
                            #low="red1", high="dodgerblue1",
                            low="gray98", high="gray30",
                            limits=c(0,  quantile(Knockoffs.window$fracs, 1)[[1]]),
                            labels = c("Quantile 0", "", "", "Quantile 1"),
                            breaks=c(0,quantile(Knockoffs.window$fracs, 0.5)[[1]],
                                     quantile(Knockoffs.window$fracs, 0.75)[[1]], 
                                     quantile(Knockoffs.window$fracs, 1)[[1]]),
                            space="Lab", na.value="gray", guide="colourbar") # ,labels=function(x){x}
      
      
    }
    
    if(!evals & !fill_in) {
      plot_chicago <- Knockoffs.window %>%
        mutate(Resolution=as.character(Resolution)) %>%
        mutate(Resolution = paste0("res", Resolution)) %>%
        mutate(Height=resolution.heights[Resolution]) %>%
        ggplot() +
        geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5), #, 
                  color="black")  
    }
    
  } else {
    plot_chicago <- ggplot(tibble()) + geom_blank()
  }
  
  plot_chicago <- plot_chicago +
    ylab("Resolution (Mb)") + xlab("") +
    coord_cartesian(xlim = c(window.left,window.right)) +
    scale_x_continuous(expand=c(0.01,0.01), labels=bp.labeler_UKB) +
    scale_y_continuous(limits=c(0.5,max(resolution.heights)+0.5),
                       labels=resolution.labels, breaks=resolution.heights) +
    #ggtitle(paste0("Chicago plot Chromosome ", window.chr)) +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          panel.border=element_blank(),
          panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
          panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray"),
          text = element_text(size=font.size),
          axis.title.y = element_text(size=title.font.size),
          plot.title = element_text(size=title.font.size),
          legend.text = element_text(size=legend.font.size),
          legend.title = element_text(size=legend.font.size),
          legend.key.height = unit(0.75,"line")
    )
  
  return(plot_chicago)
}


# PLOT FOR UKB INCLUDING IMPLICITY REJECTED RESOLUTIONS 

# Function to plot a chicago plot for UKB results including implicitly rejected groups 
# INPUT: 
  # window.chr: which chromosome to plot 
  # window.left: BP starting position for plot
  # window.right: BP end position for plot
  # Discoveries: data frame containing information on the discoveries made by a method 
  # more_info_additional_groups_implicity_rejected: data frame containing information on the implicitly rejected groups
  # evals: indicator for whether to fill based on evalues or evalues without multiplier (fractions)
# OUTPUT: 
# plot_chicago: chicago plot


plot_chicago_with_bigger_res <- function(window.chr, window.left, window.right, Discoveries,
                                         more_info_additional_groups_implicity_rejected, evals = FALSE, fill = TRUE) {
  
  # Extract knockoff discoveries within this window
  if(!is.null(Discoveries)) {
    Knockoffs.window <- Discoveries %>%
      filter(CHR==window.chr, BP.min<=window.right, BP.max>=window.left)
    #cat(sprintf("There are %d knockoff discoveries within this window.\n", nrow(Knockoffs.window)))
    
    additional.window <- more_info_additional_groups_implicity_rejected %>%
      filter(CHR==window.chr, BP.min<=window.right, BP.max>=window.left)
    
  } else {
    Knockoffs.window <- tibble()
  }
  
  # Plot knockoff discoveries
  resolution.list <- c("res0", "res1", "res2", "res3", "res4", "res5", "res6") %>% rev
  resolution.heights <- seq(length(resolution.list))
  names(resolution.heights) <- resolution.list
  #resolution.labels <- paste(parse_number(resolution.list), "\\%", sep="")
  #resolution.labels <- resolution.list
  resolution.labels <- c("single-SNP", "3 kb", "20 kb", "41 kb", "81 kb", "208 kb", "425 kb") %>% rev
  #resolution.labels <- c("single-SNP", "3", "20", "41", "81", "208", "425") %>% rev
  
  if(nrow(Knockoffs.window)>0) {
    
    plot_chicago_original_df <- Knockoffs.window %>%
      mutate(Resolution=as.character(Resolution)) %>%
      mutate(Resolution = paste0("res", Resolution)) %>%
      mutate(Height=resolution.heights[Resolution]) 
    
    
    plot_chicago_additional_df <- more_info_additional_groups_implicity_rejected %>%
      mutate(Resolution=as.character(Resolution)) %>%
      mutate(Resolution = paste0("res", Resolution)) %>%
      mutate(Height=resolution.heights[Resolution]) 
    
    
    if(fill) {
      ggplot(plot_chicago_original_df) + 
        geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5, fill=fracs - mean(fracs)),
                  color="black") +
        scale_fill_gradient(name="E-value",
                            #low="red1", high="dodgerblue1",
                            low="gray98", high="gray30",
                            limits=c(0,  quantile(Knockoffs.window$fracs, 1)[[1]]),
                            labels = c("Quantile 0", "", "", "Quantile 1"),
                            breaks=c(0,quantile(Knockoffs.window$fracs, 0.5)[[1]],
                                     quantile(Knockoffs.window$fracs, 0.75)[[1]],
                                     quantile(Knockoffs.window$fracs, 1)[[1]]),
                            space="Lab", na.value="gray", guide="colourbar") +
        geom_rect(data = plot_chicago_additional_df,  
                  aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5), fill = NA,
                  color="blue") -> plot_chicago
    } else {
      ggplot(plot_chicago_original_df) + 
        geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5),
                  color="black") +
        geom_rect(data = plot_chicago_additional_df,  
                  aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5), fill = NA,
                  color="blue") -> plot_chicago
    }
    
   
    
  } else {
    plot_chicago <- ggplot(tibble()) + geom_blank()
  }
  
  plot_chicago <- plot_chicago +
    ylab("Resolution (Mb)") + xlab("") +
    coord_cartesian(xlim = c(window.left,window.right)) +
    scale_x_continuous(expand=c(0.01,0.01), labels=bp.labeler_UKB) +
    scale_y_continuous(limits=c(0.5,max(resolution.heights)+0.5),
                       labels=resolution.labels, breaks=resolution.heights) +
    #ggtitle(paste0("Chicago plot Chromosome ", window.chr)) +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          panel.border=element_blank(),
          panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
          panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray"),
          text = element_text(size=font.size),
          axis.title.y = element_text(size=title.font.size),
          plot.title = element_text(size=title.font.size),
          legend.text = element_text(size=legend.font.size),
          legend.title = element_text(size=legend.font.size),
          legend.key.height = unit(0.75,"line")
    )
  
  return(plot_chicago)
}


# all in a single line 
plot_chicago_single_line <- function(window.chr, window.left, window.right, Discoveries, evals = FALSE, 
                                         fill_color = FALSE) {
  
  # Extract knockoff discoveries within this window
  if(!is.null(Discoveries)) {
    Knockoffs.window <- Discoveries %>%
      filter(CHR==window.chr, BP.min<=window.right, BP.max>=window.left)
    #cat(sprintf("There are %d knockoff discoveries within this window.\n", nrow(Knockoffs.window)))

  } else {
    Knockoffs.window <- tibble()
  }
  
  # Plot knockoff discoveries
  resolution.list <- c("res0", "res1", "res2", "res3", "res4", "res5", "res6") %>% rev
  resolution.heights <- seq(length(resolution.list))
  names(resolution.heights) <- resolution.list
  #resolution.labels <- paste(parse_number(resolution.list), "\\%", sep="")
  #resolution.labels <- resolution.list
  resolution.labels <- c("single-SNP", "3 kb", "20 kb", "41 kb", "81 kb", "208 kb", "425 kb") %>% rev
  #resolution.labels <- c("single-SNP", "3", "20", "41", "81", "208", "425") %>% rev
  
  if(nrow(Knockoffs.window)>0) {
    
    plot_chicago_original_df <- Knockoffs.window %>%
      mutate(Resolution=as.character(Resolution)) %>%
      mutate(Resolution = paste0("res", Resolution)) %>%
      mutate(Height=1) 
    
    
    if(!fill_color) {
      ggplot(plot_chicago_original_df) + 
        geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5),
                  color="black", fill="black")   -> plot_chicago
    } else {
      ggplot(plot_chicago_original_df) + 
        geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5, fill= fracs - mean(fracs)),
                  color="black") + 
        ggtitle("Rejected regions") +
        scale_fill_gradient(name="E-value",
                            low="gray98", high="gray30",
                            limits=c(0,  quantile(Knockoffs.window$fracs, 1)[[1]]),
                            labels = c("Quantile 0", "", "", "Quantile 1"),
                            breaks=c(0,quantile(Knockoffs.window$fracs, 0.5)[[1]],
                                     quantile(Knockoffs.window$fracs, 0.75)[[1]],
                                     quantile(Knockoffs.window$fracs, 1)[[1]]),
                            space="Lab", na.value="gray", guide="colourbar") -> plot_chicago
    }
    
   
    
   

  } else {
    plot_chicago <- ggplot(tibble()) + geom_blank()
  }
  
  plot_chicago <- plot_chicago +
    ylab("") + xlab("") +
    coord_cartesian(xlim = c(window.left,window.right)) +
    scale_x_continuous(expand=c(0.01,0.01), labels=bp.labeler_UKB) +
    scale_y_continuous(limits=c(0.5,1.5), breaks = NULL) +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          panel.border=element_blank(),
          panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
          panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray"),
          text = element_text(size=font.size),
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_text(size=title.font.size),
          plot.title = element_text(size=title.font.size),
          legend.text = element_text(size=legend.font.size),
          legend.title = element_text(size=legend.font.size),
          legend.key.height = unit(0.75,"line")
    )
  
  return(plot_chicago)
}


plot_annotations <- function(window.chr, window.left, window.right, Annotations.func) {
  
  
  # Extract color map
  annotation.color.map <- Annotations.func %>% group_by(name, itemColor) %>% summarise() %>%
    ungroup() %>%
    mutate(name.num=parse_number(as.character(name))) %>%
    mutate(label=gsub("\\d+_", "",name), label=gsub(fixed("_"), " ",label)) %>%
    mutate(label=as.factor(label)) %>%
    arrange(name.num)
  
  # Convert names to factors according to color maps
  Annotations.func <- Annotations.func %>%
    mutate(label=gsub("\\d+_", "",name), label=gsub(fixed("_"), " ",label)) %>%
    mutate(label=factor(label, levels=annotation.color.map$label, labels=annotation.color.map$label))
  
  # Plot functional annotations
  myColors <- annotation.color.map$itemColor
  names(myColors) <- annotation.color.map$label
  
  # Select functional annotations within this window
  Functional.window <- Annotations.func %>%
    filter(chrom==window.chr, chromStart<=window.right, chromEnd>=window.left)
  #cat(sprintf("There are %d functional annotations within this window.\n",
  #            nrow(Functional.window)))
  
  # Make plot
  if(nrow(Functional.window)>0) {
    p.functional <- Functional.window %>%
      mutate(chromStart=chromStart, chromEnd=chromEnd) %>%
      ggplot() +
      geom_rect(aes(xmin=chromStart, xmax=chromEnd, ymin=0.5, ymax=1.5, fill=label)) +
      ylab("") + xlab("") +
      coord_cartesian(xlim = c(window.left,window.right)) +
      scale_x_continuous(expand=c(0.01,0.01), labels=bp.labeler_UKB) +
      scale_color_manual(values=myColors, guide=FALSE) +
      scale_fill_manual(values=myColors, name="Variant annotation") +
      ggtitle("Functional annotations") +
      theme_void() +
      theme(legend.key.size = unit(0.5,"line"),
            text = element_text(size=font.size),
            plot.title = element_text(size=title.font.size),
            axis.title = element_text(size=axis.font.size),
            legend.text = element_text(size=legend.font.size),
            legend.title = element_text(size=legend.font.size),
      )+
      guides(fill=guide_legend(ncol=2))
  } else {
    p.functional <- ggplot(tibble()) + geom_blank()
  }
  
  # Return plot
  return(p.functional)
}


load_annotations <- function(data_dir) {
  # Functional annotations
  annotations.func.file <- sprintf("%s/annotations/wgEncodeBroadHmmGm12878HMM.txt", data_dir)
  Annotations.func.raw <- read_tsv(annotations.func.file, 
                                   col_names = FALSE,
                                   col_types=cols(col_integer(), 
                                                  col_character(),
                                                  col_integer(),
                                                  col_integer(),
                                                  col_character(),
                                                  col_integer(),
                                                  col_character(),
                                                  col_integer(),
                                                  col_integer(),
                                                  col_integer()))
  names(Annotations.func.raw) = c("#bin", "chrom", "chromStart", "chromEnd",
                                  "name", "score", "strand", "thickStart", 
                                  "thickEnd", "itemRgb")
  
  Annotations.func.raw <- Annotations.func.raw %>% mutate(name=as.factor(name))
  
  # Remove weird chromosomes and convert colors
  valid.chrom <- paste("chr", seq(1,22), sep="")
  Annotations.func <- Annotations.func.raw %>%
    mutate(name=as.factor(name)) %>%
    filter(chrom %in% valid.chrom) %>% mutate(chrom=parse_number(as.character(chrom))) %>%
    mutate(itemR = as.integer(floor(itemRgb/256^2) %% 256), 
           itemG = as.integer(floor(itemRgb/256) %% 256), 
           itemB = as.integer(itemRgb %% 256)) %>%
    mutate(itemColor = rgb(red=itemR, blue=itemB, green=itemG, maxColorValue=255)) %>%
    dplyr::select(c("#bin", "chrom", "chromStart", "chromEnd", "name", "score",
                    "strand", "thickStart", "thickEnd", "itemR", "itemG", "itemB", "itemColor")) %>%
    arrange(chrom, chromStart)
  
  # Extract color map
  annotation.color.map <- Annotations.func %>% group_by(name, itemColor) %>% summarise() %>% 
    ungroup() %>%
    mutate(name.num=parse_number(as.character(name))) %>% 
    mutate(label=gsub("\\d+_", "",name), label=gsub(fixed("_"), " ",label)) %>%
    arrange(name.num)
  
  # Convert names to factors according to color maps
  Annotations.func <- Annotations.func %>%
    mutate(name=factor(name, levels=annotation.color.map$name, labels=annotation.color.map$name))
  
  # Gene annotations
  annotations.genes.file <- sprintf("%s/annotations/ncbiRefSeq.txt", data_dir)
  Annotations.genes.raw <- read_tsv(annotations.genes.file, 
                                    col_names = FALSE,
                                    col_types=cols(col_integer(),
                                                   col_character(),
                                                   col_character(),
                                                   col_character(),
                                                   col_integer(),
                                                   col_integer(),
                                                   col_integer(),
                                                   col_integer(),
                                                   col_integer(),
                                                   col_character(),
                                                   col_character(),
                                                   col_integer(),
                                                   col_character(),
                                                   col_character(),
                                                   col_character(),
                                                   col_character()))
  names(Annotations.genes.raw) = c("#bin", "name", "chrom", "strand", "txStart",
                                   "txEnd", "cdsStart", "cdsEnd", "exonCount",
                                   "exonStarts", "exonEnds", "score", "name2", 
                                   "cdsStartStat", "cdsEndStat", "exonFrames")
  
  # Remove weird chromosomes
  valid.chrom <- paste("chr", seq(1,22), sep="")
  Annotations.genes <- Annotations.genes.raw %>%
    filter(chrom %in% valid.chrom) %>% mutate(chrom=parse_number(as.character(chrom))) %>%
    arrange(chrom, `#bin`)
  
  # Split rows corresponding to same gene but different exons
  Exons <- Annotations.genes %>% 
    separate_rows(exonStarts, exonEnds, exonFrames, sep=",", convert=TRUE) %>%
    drop_na()
  
  # Pick the canonical transcripts
  # That is, for each unique "name2", keep only the rows corresponding to the "name" 
  # with the largest sum of exon lengths
  Exons.canonical <- Exons %>%
    mutate(exonLength=exonEnds-exonStarts) %>%
    group_by(name, name2) %>% summarise(Length=sum(exonLength)) %>%
    ungroup() %>% group_by(name2) %>% top_n(1, Length) %>%
    inner_join(Exons, by=c("name", "name2")) 
  
  annotations = c()
  annotations$Annotations.func = Annotations.func
  annotations$Exons.canonical = Exons.canonical
  
  return(annotations)
}


plot_combined <- function(window.chr, window.left, window.right, Discoveries, 
                          data_dir,
                          highlight.gene=NULL, max.gene.rows=10) {
  
  # Make sure that the window is not empty
  if(window.right<=window.left) {
    return(ggplot(tibble()) + geom_blank())
  }
  
  annotations <- load_annotations(data_dir)
  Annotations.func <- annotations$Annotations.func
  
  # Make Chicago plot with KnockoffZoom discoveries
  p.knockoffs <- plot_chicago_single_line(window.chr, window.left, window.right, Discoveries, fill_color = TRUE)
  
  # Plot functional annotations
  if(!is.null(Annotations.func)) {
    p.functional <- plot_annotations(window.chr, window.left, window.right, Annotations.func)
  } else {
    p.functional <- ggplot(tibble()) + geom_blank()
  }
  
  # Determine relative heights of each subplot
  height.knockoffs <- 0.3
  height.functional <- 0.3
  heights <- c(height.knockoffs, height.functional)
  
  # Convert the plot objects for placement
  debug.lines <- FALSE
  g1 <-  ggplotGrob(p.knockoffs + theme(legend.position="none"))
  g2 <-  ggplotGrob(p.functional + theme(legend.position="none"))
  fg1 <- gtable_frame(g1, width = unit(1, "null"), height = unit(heights[1], "null"), debug = debug.lines)
  fg2 <- gtable_frame(g2, width = unit(1, "null"), height = unit(heights[2], "null"), debug = debug.lines)
  
  # Combine the main plots
  fg.l <- gtable_frame(gtable_rbind(fg1, fg2),
                       width = unit(4, "null"), height = unit(1, "null"))
  
  # Extract the legends
  g7 <- ggplotGrob(ggplot())
  try(g7 <- ggplotGrob(ggplotify::as.ggplot(get_legend(p.knockoffs))+
                         theme(text = element_text(size=legend.font.size))
  ), silent=FALSE)
  g8 <- ggplotGrob(ggplot())
  try(g8 <- ggplotGrob(ggplotify::as.ggplot(get_legend(p.functional))+
                         theme(text = element_text(size=legend.font.size))
  ), silent=TRUE)
  
  # Combine the legends
  g0 <- ggplotGrob(ggplot())
  fg7 <- gtable_frame(g7, width = unit(1, "null"), height = unit(heights[1], "null"), debug = debug.lines)
  fg8 <- gtable_frame(g8, width = unit(1, "null"), height = unit(heights[2], "null"), debug = debug.lines)
  fg00 <- gtable_frame(g0, width = unit(1, "null"), height = unit(1, "null"), debug = debug.lines)
  fg.r <- gtable_frame(gtable_rbind(fg7,fg8),
                       width = unit(1, "null"), height = unit(1, "null"))
  
  # Combine main plots and legends
  grid.newpage()
  combined <- gtable_frame(gtable_cbind(fg.l, fg.r),
                           width = unit(1, "null"),
                           height = unit(1, "null"))
  p.final <- ggplotify::as.ggplot(combined)
  
  # Return complete plot
  return(p.final)
}


