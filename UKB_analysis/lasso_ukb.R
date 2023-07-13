#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)

# IMPORTANT: The data is not available; we have applied for it from a source (the UKBiobank) and interested parties can also apply to the same source

# population can be one of whitenonbritish or british
# all samples are unrelated
my_population <- as.character(args[1])
pheno.name <- as.character(args[2])
resolution <- as.numeric(args[3])
ncores <- as.numeric(args[4])

print(paste0("Population: ", my_population, 
             "; Phenotype: ", pheno.name, 
             "; Resolution: ", resolution, 
             "; N Cores: ", ncores))

outdir <- "/"
dfmax  <- 10000

# Set seed for cross-validated lasso statistics
set.seed(2022)

# partially using code from https://github.com/msesia/knockoffgwas/blob/master/knockoffgwas/utils/lasso.R

# Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(bigsnpr))
suppressMessages(library(knockoff))


####################
## Load genotypes ##
####################

bed_bim_fam_filename <- paste0("INSERT_FILE_PATH")


cat("Reading genotypes ... ")
print(bed_bim_fam_filename)

plinkfile = paste0(bed_bim_fam_filename,
                   ".bed")
if(!file.exists(paste0(bed_bim_fam_filename, ".bk"))) {
  x = snp_readBed2(plinkfile)
}
cat("done.\n")


# Attach the "bigSNP" object in R session
cat("Attaching bigSNP object... ")
obj.bigSNP <- snp_attach(paste0(bed_bim_fam_filename,
                                ".rds"))
cat("done.\n")

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
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
SNP <- obj.bigSNP$map$marker.ID

#############################
## Load variant partitions ##
#############################

## Load list of variants and partitions
bim.file <- paste0(bed_bim_fam_filename, ".bim")
Variants <- read.csv(bim.file, sep = "\t", header = FALSE)
colnames(Variants) <-  c("CHR", "SNP", "X0", "BP", "X1", "X2")

Variants <- Variants %>%
  mutate(Knockoff = ifelse(endsWith(SNP, ".k"), TRUE, FALSE), 
         SNP = sub(".k$", "", SNP)) %>% 
  dplyr::select(CHR, SNP, BP, Knockoff)

# Import Group Information 
grp.file <- c()
chr <- seq(1, 22)
for(c in chr) {
  
  myfilename <- paste0("INSERT_FILE_PATH", c, "_ibd1_res", resolution, "_grp",".txt")
  myfile <- read.csv(myfilename, sep="") 
  
  myfile$CHR <- c
  
  grp.file <- rbind(grp.file, myfile)
  
}


Variants <- Variants %>% 
  mutate(CHR = as.integer(CHR)) %>%
  left_join(grp.file, by = c("SNP", "CHR"))

group_bp_min_max <- Variants %>% 
  group_by(CHR, Group) %>% 
  summarise(BP.min = min(BP), BP.max = max(BP))


# Compute scaling factor for the genotypes
cat("Computing scaling factors for all variants... ")
scaler <- big_scale()
G.scale <- scaler(G)
scaling.factors <- G.scale$scale
cat("done.\n")

#####################
## Load phenotypes ##
#####################

cat("Reading phenotype file... ")
# Load phenotype table
Phenotypes <- read.csv(paste0("INSERT_FILE_PATH", 
                              my_population, "only" ,".reordered.csv"))
cat("done.\n")



y <- Phenotypes[, pheno.name]


# Find the class of response (numeric or binary factor)
y.unique <- unique(y[!is.na(y)])
if(length(y.unique)==2) {
  phenotype.class <- "binary"
  y <- factor(y, levels=c(0,1), labels=c(0,1))
  y <- as.numeric(levels(y))[y]
} else {
  phenotype.class <- "continuous"
}

###################
## Fit the lasso ##
###################

# extract covariates 
Covariates <- Phenotypes %>% dplyr::select(age, age_squared, sex, PC1, PC2, PC3, PC4, PC5)
Covariates <- as.matrix(Covariates)

cat("Running lasso... ")
# Fit the lasso
if(phenotype.class=="binary") {
  cat(sprintf("Fitting sparse logistic regression with %d observations, %d variants and %d covariates... ",
              length(y), ncol(G), ncol(Covariates)))
  lasso.fit <- big_spLogReg(G, 
                            y01.train = y, 
                            covar.train=Covariates,
                            pf.covar = rep(0, ncol(Covariates)), # don't penalize covariates
                            dfmax=dfmax, 
                            ncores=ncores)
} else {
  cat(sprintf("Fitting sparse linear regression with %d observations, %d variants and %d covariates... ",
              length(y), ncol(G), ncol(Covariates)))
  lasso.fit <- big_spLinReg(G, 
                            y.train=y,
                            covar.train=Covariates,
                            pf.covar = rep(0, ncol(Covariates)), # don't penalize covariates
                            dfmax=dfmax, 
                            ncores=ncores)
}
cat("done.\n")


# 10-fold CV (default)
# nabort = 10 (default)
# nlam.min = 50 (default)

# Extract beta from each fold and combine them
cat("Extracting regression coefficients... ")
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
  select(CHR, SNP, BP, Knockoff, Z)

# Extract the estimated coefficients
Lasso.res <- Beta %>%
  inner_join(Variants, by = c("CHR", "SNP", "BP", "Knockoff")) %>%
  filter(Z!=0) %>%
  arrange(desc(Z)) %>%
  select(CHR, SNP, BP, Z, Group, Knockoff)
cat("done.\n")


##################################
## Compute the test statistics ##
##################################

cat("Computing test statistics... ")

# Compute the knockoff statistics
W.stats <- function(Z, knockoff) {
  importance <- abs(Z)
  z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
  zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
  w <- z-zk
}


Stats <- Lasso.res %>%
  select("CHR", "Group", "SNP", "BP", "Knockoff", "Z") %>%
  filter(Z!=0) %>%
  group_by(CHR, Group) %>%
  summarize(W = W.stats(abs(Z),Knockoff),
            Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
            Size=n()) %>%
  ungroup() %>%
  arrange(desc(abs(W))) %>%
  select(CHR, Group, SNP.lead, BP.lead, Size, W) %>%
  filter(W!=0)

# add in BP min and BP max information
Stats <- left_join(Stats, group_bp_min_max, by = c("CHR", "Group"))

cat("done.\n")

# Save results
stats.file <- paste0(outdir, "/" , "lasso_",pheno.name,"_",
                     my_population,  "_unrelated_res", resolution, ".txt")
Stats %>% write_delim(stats.file, delim=" ")
cat(sprintf("Test statistics written on:\n%s\n", stats.file))

saveRDS(lasso.fit, file=paste0(outdir, "/" , pheno.name,"_",
                               my_population,  "_res", resolution,"_lasso_fit"))
write.csv(Lasso.res, file=paste0(outdir, "/" , pheno.name,"_", 
                                 my_population,  "_res", resolution, "_lasso_beta"))
