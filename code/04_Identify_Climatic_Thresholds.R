# Identify_Climatic_Thresholds.R
# ---------------------------------------------------------------------------------- 
# THIS SCRIPT ESTIMATES AND RECORDS THE BREAKPOINT IN THE NONLINEAR RELATIONSHIPS  
# BETWEEN CLIMATE AND THE 30-YR PROBABILITY OF FIRE OCCURENCE. R version 3.2.2 
# (2015-08-14) WAS USED TO CONDUCT THIS ANALYSIS.                                       
#                                                                                   
# FILE REQUIREMENTS:
#   (1) Partial dependence plots for Temperature of the warmest month (TempWarm) and 
#       annual moisture availability (AnnDEF). In each model results directory (e.g.,
#       'AK_FINAL_RESULTS') in the "Generalized Boosted Models and Analysis Scripts." 
#       These files are named:
#
#       - 'partDep_TempWarm.csv'
#       - 'partDep_AnnDEF.csv'
#
# DEPENDENCIES:
#   * segmented package (version 0.5-1.4)
#
# CITATION:
# Young, A. M. et al. Climatic thresholds shape northern high-latitude fire 
# regimes and imply vulnerability to future climate change. 2017.
# DOI:10.1111/ecog.02205
#
# Author(s): Adam Young                                                             
# Date created: March 2014                                                              
# Edited for publication: January 2016   
#
# Contact: Philip Higuera, PhD, philip.higuera[at]umontana.edu                                                    
# ---------------------------------------------------------------------------------- 

# INITIALIZE WORKSPACE
rm(list = ls()) # CLEAR WORKSPACE 
graphics.off() # DELETE CURRENT GRAPHICS WINDOWS
cat("\14") # CLEAR COMMAND PROMPT

wdir <- "/.../" # Parent (working) directory where all subfolders from 
                # nested datasets are stored. This may need to be further 
                # modified below depending on how directories are organized 
                # by user.

mdl_dir <- c("AK_FINAL_RESULTS","BOREAL_FINAL_RESULTS","TUNDRA_FINAL_RESULTS")

# LOAD REQUIRED LIBRARIES 
require(segmented)

# THE FOLLOWING LIST PROVIDES MATRICES INITIAL ESTIMATES WHERE THE BREAKPOINTS IN THE
# NONLINEAR RELATIONSHIP OCCUR, AS WELL AS LIMITS TO THE RANGE OF COVARIATE VALUES 
# SHOULD BE LIMITED TO. THESE VALUES WERE OBTAINED FROM VISUAL INSPECTION OF RESULTS
# FROM BRT ANALYSIS.
bps_twrm <- list(ak = matrix(c(14,12,17),1,3),
                 boreal = matrix(c(14,12,17),1,3),
                 tundra = matrix(c(13,8,15),1,3))

bps_adef <- list(ak = matrix(c(200,-100,300),1,3),
                 boreal = matrix(c(200,0,300),1,3),
                 tundra = matrix(c(-200,-300,-50,200,-50,300),2,3,byrow = TRUE))

nboot <- 2000 # NUMBER OF BOOTSTRAP SAMPLES TO RECORD

for (m in 1:length(mdl_dir)){ # FOR EACH MODEL ...
  
  setwd(paste(wdir,mdl_dir[m],sep = "")) # ... CHANGE THE WORKING DIRECTORY
  
  # TEMPWARM CLIMATIC THRESHOLDS
  partDep <- read.csv("partDep_TempWarm.csv") # LOAD RESULTS FROM PARTIAL DEPENDENCE
  bps_i <- bps_twrm[[m]] # INITIAL ESTIMATED BREAKPOINT LOCATIONS
  nbps <- nrow(bps_i) # NUMBER OF BREAKPOINTS TO ESTIMATE
  boot_bp <- matrix(NA,nboot,nbps) # CREATE MATRIX TO STORE RESULTS
  for (n in 1:nbps){
    bp_i <- bps_i[n,1] # CURRENT ESTIMATED BREAKPOINT TO ESTIMATE
    variable <- partDep[,2] # COVARIATE VALUES
    # LIMIT COVARIATE VALUES TO THOSE IMMEDIATELY AROUND THE ESTIMATED BREAKPOINT
    variable[(variable < bps_i[n,2] | variable > bps_i[n,3])] <- NA
    for (i in 1:nboot){ # FOR EACH BOOTSTRAP SAMPLE ...
      set.seed(i)
      # RANDOMLY SAMPLE THE PARTIAL DEPENDENCE RESULTS FROM THE 100 BRTS WITH
      # REPLACEMENT.
      samps <- sample(3:ncol(partDep),size = length(3:ncol(partDep)),replace = TRUE)
      # CALUCATE THE MEDIAN PREDICTED PROBABILTIY OF FIRE OCCURRENCE FOR THE 
      # COVARIATE VALUES OF INTEREST FROM THE 100 RANDOMLY SAMPLED PARTIAL 
      # DEPENDENCE RESULTS
      predprob_i <- apply(partDep[,samps],MARGIN = 1,FUN = median)
      # SET UP TEMPORARY DATA FRAME TO STORE THESE DATA - FOR USE IN SEGMENTED
      # REGRESSION
      dataframe_i <- data.frame(y = predprob_i,x = variable)
      # CONDUCT SEGMENTED REGRESSION ANALYSIS
      seg_lm <- lm(y ~ x,data = dataframe_i)
      est_bp <- segmented(seg_lm,seg.Z = ~ x,psi = bp_i)
      # RECORD BREAKPOINT BOOTSTRAP ESTIMATE
      boot_bp[i,n]   <- est_bp$psi[2]
    }
  }
  write.csv(boot_bp,"climThresholds_TempWarm.csv") # EXPORT RESULTS
  
  # ANNDEF CLIMATIC THRESHOLDS - SAME AS THOSE ABOVE FOR 
  partDep <- read.csv("partDep_AnnDEF.csv")
  bps_i <- bps_adef[[m]] # INITIAL ESTIMATED BREAKPOINT LOCATIONS
  nbps <- nrow(bps_i) # NUMBER OF BREAKPOINTS TO ESTIMATE
  boot_bp <- matrix(NA,nboot,nbps) # CREATE MATRIX TO STORE RESULTS
  for (n in 1:nbps){
    bp_i <- bps_i[n,1] # CURRENT ESTIMATED BREAKPOINT TO ESTIMATE
    variable <- partDep[,2] # COVARIATE VALUES
    # LIMIT COVARIATE VALUES TO THOSE IMMEDIATELY AROUND THE ESTIMATED BREAKPOINT
    variable[(variable < bps_i[n,2] | variable > bps_i[n,3])] <- NA
    for (i in 1:nboot){ # FOR EACH BOOTSTRAP SAMPLE ...
      set.seed(i)
      # RANDOMLY SAMPLE THE PARTIAL DEPENDENCE RESULTS FROM THE 100 BRTS WITH
      # REPLACEMENT.
      samps <- sample(3:ncol(partDep),size = length(3:ncol(partDep)),replace = TRUE)
      # CALUCATE THE MEDIAN PREDICTED PROBABILTIY OF FIRE OCCURRENCE FOR THE 
      # COVARIATE VALUES OF INTEREST FROM THE 100 RANDOMLY SAMPLED PARTIAL 
      # DEPENDENCE RESULTS
      predprob_i <- apply(partDep[,samps],MARGIN = 1,FUN = median)
      # SET UP TEMPORARY DATA FRAME TO STORE THESE DATA - FOR USE IN SEGMENTED
      # REGRESSION
      dataframe_i <- data.frame(y = predprob_i,x = variable)
      # CONDUCT SEGMENTED REGRESSION ANALYSIS
      seg_lm <- lm(y ~ x,data = dataframe_i)
      est_bp <- segmented(seg_lm,seg.Z = ~ x,psi = bp_i)
      # RECORD BREAKPOINT BOOTSTRAP ESTIMATE
      boot_bp[i,n]   <- est_bp$psi[2]
    }
  }
  write.csv(boot_bp,"climThresholds_AnnDEF.csv") # EXPORT RESULTS
}
