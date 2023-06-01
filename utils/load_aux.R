
######################################################################

              # Load Packages and Functions #

######################################################################

# load packages
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(gglasso))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(knockoff))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(JuliaCall))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(stargazer))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dirmult))

# load functions 
source(paste0(mydir, "utils/ebh_aux.R"))
source(paste0(mydir, "utils/kelp_aux.R"))
source(paste0(mydir, "utils/filter_aux.R"))
source(paste0(mydir, "utils/knockoff_aux.R"))
source(paste0(mydir, "utils/multiplier_tuning_aux.R"))
source(paste0(mydir, "utils/plot_aux.R"))
source(paste0(mydir, "utils/power_fdr_aux.R"))
source(paste0(mydir, "utils/read_process_UKB_lasso_aux.R"))
source(paste0(mydir, "utils/sim_data_aux.R"))
source(paste0(mydir, "utils/tune_gamma_sim_aux.R"))
source(paste0(mydir, "utils/emlkf_aux.R"))

















