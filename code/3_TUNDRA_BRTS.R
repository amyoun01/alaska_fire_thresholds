# TUNDRA_BRTS.R
# ---------------------------------------------------------------------------------- 
# THIS SCRIPT USES HISTORICAL FIRE, CLIMATE, VEGETATION, AND TOPOGRAPHIC DATA TO    
# CONSTRUCT BOOSTED REGRESSION TREE MODELS, CAPTURING HISTORICAL FIRE-CLIMATE       
# RELATIONSHIPS IN ALASKA TUNDRA FROM 1950-2009. THE SCRIPT USES THE gbm, raster, and      
# R.utils PACKAGES TO COMPLETE THIS ANALYSIS. R version 3.2.2 (2015-08-14) WAS USED
# TO CONDUCT THIS ANALYSIS.
#                                                                                   
# FILE REQUIREMENTS:
#
#   Available in the "Vegetation classification and topography of Alaska" nested 
#   dataset. 
#   ---------------------------------------------------------------------------------
#       Spatial data to train and test BRT analysis:
#       - 'AK_VEG.tif': Vegetation classification per pixel
#       - 'TR.tif': Topographic ruggedness
#       - 'ecor.tif': Ecoregion classification
#   ---------------------------------------------------------------------------------
#
#   Available in the "/AncillaryData/train_test_data/" folder within the "Generalized 
#   Boosted Models and Analysis Scripts" nested dataset. 
#   ---------------------------------------------------------------------------------
#       - Climate normals constructed from randomized 30-yrs into training and
#         testing datasets.
#       - Fire frequency maps for randomized thirty-year periods to test and train
#         BRTs.
#   ---------------------------------------------------------------------------------
#
# DEPENDENCIES:
#   * aucroc.R - CUSTOM WRITTEN FUNCTION THAT DOES THE AUC ANALYSIS
#   * gbm package (version 2.1.1)
#   * raster package (version 2.4-20)
#   * R.utils package (version 2.2.0)
#
# CITATION:
# Young, A. M. et al. Climatic thresholds shape northern high-latitude fire 
# regimes and imply vulnerability to future climate change. 2017.
# DOI:10.1111/ecog.02205
#
# Author(s): Adam Young                                                             
# Date created: May 2013                                                              
# Edited for publication: January 2016   
#
# Contact: Philip Higuera, PhD, philip.higuera[at]umontana.edu                                                    
# ---------------------------------------------------------------------------------- 

#-----------------------------------------------------------------------------------#
########################### INITIALIZE WORKSPACE ####################################
#-----------------------------------------------------------------------------------#

rm(list = ls()) # CLEAR WORKSPACE 
graphics.off() # DELETE CURRENT GRAPHICS WINDOWS
cat("\14") # CLEAR COMMAND PROMPT

wdir <- "/.../" # Parent (working) directory where all subfolders from 
                # nested datasets are stored. This may need to be further 
                # modified below depending on how directories are organized 
                # by user.

save_dir <- paste(wdir,"TUNDRA_FINAL_RESULTS/",sep="")

# CREATE DIRECTORY, IF IT DOES NOT ALREADY EXIST, TO STORE RESULTS
direxist <- dir.exists(save_dir)
if (direxist == FALSE){
  dir.create(save_dir,
             recursive = TRUE)
}

#-----------------------------------------------------------------------------------#
################# LOAD REQUIRED LIBRARIES AND DATASETS ##############################
#-----------------------------------------------------------------------------------#

require(gbm) # gbm: GENERALIZED BOOSTING MODELS PACKAGE
require(raster) # raster PACKAGE
require(R.utils) # R.utils PACKAGE

# LOAD CUSTOM FUNCTION(S)
source(paste(wdir,"Scripts/Functions/aucroc.R",
             sep = ""))# FUNCTION FOR AUC CALCULATION

#LOAD VEGETATION, ECOREGION, AND TOPOGRAPHIC SPATIAL DATA AS VECTORS USING THE
# raster() AND getValues() FUNCTIONS FROM THE raster PACKAGE. WILL HAVE TO MODIFY THE 
# DIRECTORY NAME TO MATCH THE LOCATION OF THESE DATA.
setwd(paste(wdir,"VegLandscapeData/",
            sep = ""))
mapInfo <- raster("AK_VEG.tif") # GEOGRAPHIC INFORMATION NEEDED FOR EXPORTING GEOTIFF
                                # FILES
akVeg   <- getValues(mapInfo) # [categorical, 1-5] EXPLANATORY
TR      <- getValues(raster("TR.tif")) # [meters] EXPLANATORY
ecor    <- getValues(raster("ecor.tif")) # [categorical, 1-22]

# REMOVE BOREAL ECOREGIONS FROM ANALYSIS
# LOAD ECOREGOIN CLASSIFICATION INFORMATION
ecorinfo <- read.csv(paste(wdir,"AncillaryData/ecor_info.csv",
                           sep = ""),
                     header = FALSE)
for (i in 1:nrow(ecorinfo)){
 if (ecorinfo[i,3] == 2){
  ecor[ecor == i] <- NA 
 }
}

# CLASSIFY PIXELS NOT IN STUDY AREA AS 'NA' VALUES
ecor[ecor == -9999] <- NA
akVeg[(akVeg == -9999 | is.na(ecor) == TRUE | akVeg == 5)] <- NA
ecor[is.na(akVeg) == TRUE] <- NA

# CREATE A VECTOR OF UNIQUE ECOREGION VALUES [1 - 22]
ecor_idx <- sort(unique(ecor)) # OBTAIN VECTOR OF UNIQUE ECOREGIONS INDICES

#-----------------------------------------------------------------------------------#
################### INITIALIZE AND DECLARE VALUES FOR ANALYSIS ######################
#-----------------------------------------------------------------------------------#

nmodels <- 100 # NUMBER OF MODELS TO RUN
nobs    <- 5526 # NUMBER OF PIXELS TO RANDOMLY SELECT FROM STUDY DOMAIN TO CONDUCT
                # ANALYSIS

# META-PARAMETERS FOR RUNNING BOOSTED REGRESSION TREE MODELS
n_trees           <- 5000 # NUMBER OF TREES
shrinkage         <- 0.01 # LEARNING RATE (SAME AS SHRINKAGE)
interaction_depth <- 2 # NUMBER OF PAIRWISE INTERACTIONS PER TREE
bag_fraction      <- 0.5 # BAGGING FRACTION
train_frac        <- 1.0 # TRAINING FRACTION
n_minobsinnode    <- 1 # MINIMUM NUMBER OF OBSERVATIONS ALLOWED IN EACH TREE NODE 
cv_folds          <- 5 # NUMBER OF CROSS VALIDATION PARTITIONS
verbose           <- FALSE # REPORT BOOSTED REGRESSION TREE PROGRESS TO COMMAND PROMPT?

# INFORMATION TO RECORD AND STORE PARTIAL DEPENDENCE RESULTS
res   <- 100 # RESOLUTION OF CONTINUOUS COVARIATES FOR PARTIAL DEPENDENCE PLOTS
twrm  <- c(2.4,20.0) # RANGE OF TWARM FOR PARTIAL DEPENDENCE PLOTS
adef  <- c(-492.4,750.0) # RANGE OF P-PET_ANN FOR PARTIAL DEPENDENCE PLOTS
tr    <- c(0,140) # RANGE OF TR FOR PARTIAL DEPENDENCE PLOTS
n_var <- 4 # NUMBER OF COVARIATES USED TO RUN BRT MODELS
n_veg <- length(unique(akVeg))-1 # NUMBER OF LEVELS IN VEGETATION FACTOR

#-----------------------------------------------------------------------------------#
######### CREATE PRE-ALLOCATED SPACE FOR RELEVENT DATA TO OUTPUT AS CSV FILES #######
#-----------------------------------------------------------------------------------#

# SPACE TO STORE FRP VALIDATION RESULTS
frp_t   <- matrix(NA,nrow = length(ecor_idx)*nmodels,ncol = 2)
frp_o   <- matrix(NA,nrow = length(ecor_idx),ncol = nmodels,1)
frp_p   <- matrix(NA,nrow = length(ecor_idx),ncol = nmodels,1)
# AUC AND CLASSIFICATION ACCURACY
AUC_i           <- matrix(NA,nmodels,1)
class_rates     <- matrix(NA,nrow = nmodels,ncol = 4)
thresh          <- matrix(NA,nrow = nmodels,ncol = 1)
# RELATIVE INFLUENCE OF EXPLANATORY VARIABLES RESULTS
relInf          <- matrix(NA,nrow = nmodels,ncol = n_var)
# BRT OPTIMIZATION RESULTS
nTrees          <- matrix(NA,nrow = nmodels,ncol = 1)
model_train_err <- matrix(NA,nrow = n_trees,ncol = nmodels)
model_valid_err <- matrix(NA,nrow = n_trees,ncol = nmodels)
# PARTIAL DEPENDENCE RESULTS
partDep_1  <- matrix(NA,nrow = res,ncol = nmodels+1)
partDep_2  <- matrix(NA,nrow = res,ncol = nmodels+1)
partDep_3  <- matrix(NA,nrow = res,ncol = nmodels+1)
partDep_4  <- matrix(NA,nrow = n_veg,ncol = nmodels+1)
TempWarm_AnnDEF_int  <- matrix(NA,res^2,nmodels+2)

#-----------------------------------------------------------------------------------#
################### RUN BOOSTED REGRESSION TREE ANALYSIS ############################
#-----------------------------------------------------------------------------------#

for (i in 1:nmodels){
  
  # SET SEED FOR REPEATABILITY
  set.seed(i)
  # CHANGE DIRECTORY TO WHERE TRAINING/TESTING DATA ARE STORED
  setwd(paste(wdir,"AncillaryData/train_test_data",sep = ""))
  
  # LOAD TRAINING DATA
  # LOAD FIRE FREQUENCY DATA
  fire <- getValues(raster(paste("train_firefreq_",i,".tif",sep = "")))
  fire[fire > 0] <- 1 # IF FIRE FREQUENCY IS GRETATER THAN 0, SET EQUAL TO 1. THIS WILL
                      # CREATE A PRESENCE ABSENCE MAP OF FIRE IN ALASKA.
  # MEAN TEMPERATURE OF THE WARMEST MONTH
  TempWarm <- getValues(raster(paste("train_TempWarm_",i,".tif",sep = ""))) # [degrees C]
  TempWarm[TempWarm == -9999] <- NA # SET MISSING VALUES AS 'NA'
  # TOTAL ANNUAL MOISTURE AVAILABILITY
  AnnDEF <- getValues(raster(paste("train_AnnDEF_",i,".tif",sep = ""))) # [mm]
  AnnDEF[AnnDEF == -9999] <- NA # SET MISSING VALUES AS 'NA'
  
  # CREATE DATA FRAME FOR TRAINING DATA
  train_data <- data.frame(fire     = fire,
                           TempWarm = TempWarm,
                           AnnDEF   = AnnDEF,
                           TR       = TR,
                           veg      = factor(akVeg))
  
  # LOAD TEST DATA
  # LOAD FIRE FREQUENCY DATA. HERE WE SET TWO VARIABLES (firefreq AND fire) AS WE
  # NEED THE firefreq VARIABLE TO CONSTRUCT OBSERVED FIRE ROTATION PERIODS (FRPS)
  # TO COMPARE WITH PREDICTED FRPS. THIS IS ONLY NEEDED FOR THE TESTING SET AND NOT
  # THE TRAINING SET
  firefreq <- fire <- getValues(raster(paste("test_firefreq_",i,".tif",sep = ""))) 
  fire[fire > 0] <- 1 # IF FIRE FREQUENCY IS GRETATER THAN 0, SET EQUAL TO 1. THIS 
                      # WILL CREATE A PRESENCE/ABSENCE MAP OF FIRE IN ALASKA.
  # MEAN TEMPERATURE OF THE WARMEST MONTH
  TempWarm <- getValues(raster(paste("test_TempWarm_",i,".tif",sep = ""))) # [degrees C]
  TempWarm[TempWarm == -9999] <- NA # SET MISSING VALUES AS 'NA'
  # TOTAL ANNUAL MOISTURE AVAILABILITY
  AnnDEF   <- getValues(raster(paste("test_AnnDEF_",i,".tif",sep = ""))) # [mm]
  AnnDEF[AnnDEF == -9999] <- NA # SET MISSING VALUES AS 'NA'
  
  # INITIALIZE TESTING DATA FRAME
  all_test_data  <- data.frame(fire     = fire,
                               TempWarm = TempWarm,
                               AnnDEF   = AnnDEF,
                               TR       = TR,
                               veg      = factor(akVeg))
  test_data <- all_test_data
  
  # COMPLETE CASES AND REMOVE 'NA' OBSERVATIONS FROM DATA FRAMES
  train_data       <- train_data[complete.cases(train_data$veg), ]
  test_data        <- test_data[complete.cases(test_data$veg),  ]
  
  # RANDOMLY SELECT PIXELS (POINTS) FOR TRAINING AND TESTING BRTS
  pts_train <- sample(x = 1:nrow(train_data),size = nobs,replace = FALSE)
  pts_test  <- sample(x = 1:nrow(train_data),size = nobs,replace = FALSE)
  
  # CREATE DATA FRAMES FOR TRAINING AND TESTING DATASETS THAT WILL BE USED TO CREATE 
  # AND TEST BRT ANALYSIS. TESTING DATASET WILL BE UTILIZED IN AN ROC ANALYSIS
  # (E.G., AUC).
  FINAL_train_data <- train_data[pts_train,]
  FINAL_test_data <- test_data[pts_test,]
  
  # CONSTRUCT BOOSTED REGRESSION TREE MODEL FROM TRAINING DATASET
  # THIS USES gbm VERSTION 2.1.1 RELEASED ON 3-11-2015. META-PARAMETER VALUES FOR 
  # RUNNING BRT CAN BE FOUND BEFORE FOR LOOP IN INITIALIZATION SECTION OF SCRIPT.
  brt_i <- gbm(fire ~ TempWarm + AnnDEF + TR + veg, # MODEL FORMULA
               data              = FINAL_train_data, # DATASET TO USE
               var.monotone      = NULL, # NO MONONTONIC RESTRICTIONS
               distribution      = "bernoulli", # BERNOULLI DISTRIBUTION FOR RESPONSE
               n.trees           = n_trees, # MAXIMUM NUMBER OF TREES TO TRAIN MODEL
               shrinkage         = shrinkage, # REGULARIZATION PARAMETER
               interaction.depth = interaction_depth, # TWO-WAY INTERACTIONS
               bag.fraction      = bag_fraction, # SUBSAMPLING FRACTION IN EACH ITERATION
               train.frac        = train_frac, # FRACTION OF SUPPLIED DATA TO USE 
               n.minobsinnode    = n_minobsinnode, # MINIMUM NUMBER OF OBSERVATIONS
               cv.folds          = cv_folds, # USE 5-FOLD CROSS VALIDATION
               verbose           = verbose) # PRINT PROGRESS TO SCREEN? [YES/NO]
  
  # IDENTIFY OPTIMAL NUMBER OF REGRESSION TREES (ITERATIONS) FOR iTH BRT USING THE
  # CV VALIDATION METHOD
  best_iter <- gbm.perf(brt_i,
                        method  = "cv",
                        plot.it = verbose)
  
  # CREATE A PREDICTED PROBABILTY OF FIRE MAP USING THE ALL PIXELS AVAILABLE IN THE 
  # TESTING DATA FRAME, NOT THE SUBSAMPLED DATA FRAME (I.E., FINAL_test_data).
  # WE USE ALL AVAILABLE DATA SO WE CAN CREATE MAPS FOR THE ENTIRE SPATIAL DOMAIN AND 
  # CALCULATE PREDICTED FRPS. FINAL_test_data IS USED TO EVALUATE AUC AND CLASSIFICATION
  # RATES.
  pred_map  <- predict(object  = brt_i,
                       newdata = all_test_data,
                       n.trees = best_iter,
                       type    = "response") # RETURN PROBABILITY VALUES, NOT VALUES
                                             # ON THE SCALE OF THE LINK FUNCTION
                                             # (HERE LINK = "LOGIT")
  
  # CONVERT PREDICTED PROBABILITY VALUES TO A RASTER DATA FORMAT FOR THE STATE OF AK. 
  # THIS WILL BE EXPORTED AS A .TIF FILE.
  pred_map  <- raster(t(matrix(pred_map,ncol(mapInfo),nrow(mapInfo))),
                      xmn = xmin(mapInfo),
                      xmx = xmax(mapInfo),
                      ymn = ymin(mapInfo),
                      ymx = ymax(mapInfo),
                      crs = mapInfo@crs)
  
  # THESE ARE THE PREDICTED VALUES FOR JUST THE SUBSAMPLED TESTING DATA. THIS ALLOWS
  # US TO EVALUATE AUC AND CLASSIFICATION RATES FOR THE TESTING DATASET WHILE ONLY
  # USING THE SAME NUMBER OF OBSERVATIONS USED TO TRAIN THE BRTS.
  pred  <- predict(object  = brt_i,
                   newdata = FINAL_test_data,
                   n.trees = best_iter,
                   type    = "response") # RETURN PROBABILITY VALUES, NOT VALUES
                                         # ON THE SCALE OF THE LINK FUNCTION
                                         # (HERE LINK = "LOGIT")
  
  obs   <- FINAL_test_data$fire # OBSERVED FIRE OCCURRENCE IN TESTING DATA
  
  obs_pred  <- cbind(obs,pred) # PUT OBSERVED AND PREDICTED PROBABILITIES IN SAME
                               # TWO-COLUMN MATRIX.
  # COMPLETE CASES WITH BOTH COLUMNS (I.E., REMOVE ROWS THAT HAVE MISSING VALUES)
  obs_pred  <- obs_pred[complete.cases(obs_pred[,1]), ]
  obs_pred  <- obs_pred[complete.cases(obs_pred[,2]), ]
  # FIND THOSE PREDICTED PROBABILITY VALUES ASSOCIATED WITH OBSERVED FIRE PRESENCE
  pres      <- obs_pred[obs_pred[,1] == 1, 2]
  # FIND THOSE PREDICTED PROBABILITY VALUES ASSOCIATED WITH OBSERVED FIRE ABSENCE
  abs       <- obs_pred[obs_pred[,1] == 0, 2]
  # RUN ROC ANALYSIS WICH WILL RETURN CLASSIFICATION RATES AND AUC RESULT
  roc_results = aucroc(pres, abs, n = 500)
  
  # DETERMINE THRESHOLDS AND CLASSIFCATION RATES. FIRST CALCULATE A PROBABILTY 
  # THRESHOLD THAT BEST DISCRIMINATES BETWEEN PRESENCES AND ABSENCES. HERE THAT IS
  # SUM OF THE TPR AND TNR THAT IS MAXIMIZED AT A CERTAIN PROBABILTY VALUE. max_sep
  # BELOW IS AN INDEX OF THIS PROBABILITY VALUE.
  max_sep <- which((roc_results$TPR+roc_results$TNR) == 
                    max(roc_results$TPR+roc_results$TNR))
  thr <- roc_results$thresh[max_sep] # FIND INDEXED PROBABILITY VALUE THAT BEST
                                     # DISCRIMINATES BETWEEN PRESENCE AND ABSENCES.
  # IF THE THRESHOLD IS REPRESENTED BY MULTIPLE VALUES (OR EQUAL AMONG A FEW DIFFERENT
  # PROBABILITY VALUES) TAKE THE MEAN OF SEPARATOR INDICES (I.E., max_sep).
  if (length(thr)>1){
    max_sep <- floor(mean(max_sep))
    thr <- mean(thr,na.rm = TRUE)  
  }
  
  # FIND THE VALUES OF THE CONFUSION MATRIX ASSOCATED WITH THIS PROBABILTIY THRESHOLD
  # THAT SEPARATES PRESENCES AND ABSENCES.
  class_rates[i,1] <- roc_results$TPR[max_sep] # TRUE POSITIVE RATE
  class_rates[i,2] <- roc_results$FPR[max_sep] # FALSE POSITIVE RATE
  class_rates[i,3] <- roc_results$FNR[max_sep] # FALSE NEGATIVE RATE
  class_rates[i,4] <- roc_results$TNR[max_sep] # TRUE NEGATIVE RATE
  thresh[i,1]      <- thr # IDENTIFIED PROBABILITIY THRESHOLD THAT MAXIMIZES THE
                          # SUMMATION OF THE TPR AND TNR
  
  # RECORD AUC VALUES
  AUC_i[i,1] <- roc_results$AUC
  
  # CALCULATE AND RECORD PREDICTED AND OBSERVED FRPS 
  pred_vals   <- getValues(pred_map) # PREDICTED PROBABILITY VALUES
  frp         <- matrix(NA,nrow = length(ecor_idx),ncol = 2) # STORAGE ALLOCATION FOR 
                                                             # PREDICTED AND
                                                             # OBSERVED FRP VALUES
  
  # FOR EACH ECOREGION IN STUDY DOMAIN FIND THE TOTAL NUMBER OF PIXELS THAT DEFINE
  # THE AREA OF INTEREST (I.E., EXCLUDE MISSING VALUES). THEN USE FIRE FREQUENCY
  # (I.E., INCLUDE PIXELS THAT EXPERIENCE MORE THAN ONE FIRE DURING THE TEMPORAL
  # PERIOD OF INTEREST) TO CALCULATE OBSERVED FRP. PREDICTED PROBABILTY VALUES
  # ARE USED TO CALCULATE THE PREDICTED FRP.
  for (a in 1:length(ecor_idx)){
    idx_e    <- which(ecor == ecor_idx[a])
    frp[a,1] <- 30/(sum(firefreq[idx_e],na.rm = TRUE)/length(idx_e))
    frp[a,2] <- 30/(sum(pred_vals[idx_e],na.rm = TRUE)/length(idx_e))      
  }
  
  # IF NO BURNING IS RECORED IN AN ECOREGION THEN COUNT THE OBSERVATION AS MISSING
  frp[is.infinite(frp) == TRUE] <- -9999
  # RECORD PREDICTED AND OBSERVED FRP VALUES IN TWO DIFFERENT MATRIX FORMATS.
  # THE FIRST FORMAT RECORDS ALL OBSERBED AND PREDICTED IN JUST TWO-COLUMNS. THIS IS
  # DONE FOR EASE OF PLOTTNG. THE OTHER FORMAT IS TWO SEPARATE MATRICES
  # (ONE FOR PREDICTED AND ONE FOR OBSERVED) WHERE EACH ROW REPRESENTS VALUES FOR A 
  # SPECIFIC ECOREGION. THIS IS DONE TO MORE EASILY CALCULATE ECOREGION SPECIFIC
  # FRP VALUES AND STATISTICS (E.G., THE MEDIAN OR MEAN).
  frp_t[(((i-1)*length(ecor_idx))+1):((i)*length(ecor_idx)),] <- frp
  frp_o[,i] <- frp[,1]
  frp_p[,i] <- frp[,2]
  
  # RECORD GBM OPTIMIZATION INFORMATION
  nTrees[i,1] <- best_iter
  model_train_err[,i] <- brt_i$train.error
  model_valid_err[,i] <- brt_i$cv.error
  
  # RECORD RELATIVE INFLUENCE OF COVARIATES
  varinf <- summary(brt_i,
                    n.trees = best_iter,
                    plotit = FALSE)
  I <- order(varinf$var)
  relInf[i,] <- varinf$rel.inf[I]
  
  # PARTIAL DEPENDENCE PLOT DATA
  # THESE LINES CHANGE THE VARIABLE MIN AND MAX INFORMATION (I.E., VAR.LEVELS). THIS
  # ALLOWS FOR CONSISITENT EXPORT OF PARTIAL DEPENDENCE RESULTS FOR EACH OF THE 
  brt_i$var.levels[[1]][1]  <- twrm[1]
  brt_i$var.levels[[2]][1]  <- adef[1]
  brt_i$var.levels[[3]][1]  <- tr[1]
  brt_i$var.levels[[1]][length(brt_i$var.levels[[1]])] <- twrm[2]
  brt_i$var.levels[[2]][length(brt_i$var.levels[[2]])] <- adef[2]
  brt_i$var.levels[[3]][length(brt_i$var.levels[[3]])] <- tr[2]
  
  # IF IT IS THE FIRST ITERATION OF THE FOR LOOP FIRST RECORD THE VALUES OF THE 
  # COVARIATES IN THE FIRST COLUMN OF THE MATRICES ALLOCATED TO STORING THE 
  # PARTIAL DEPENDENCE RESULTS.
  if (i==1){
    partDep_1[,1] <- plot.gbm(brt_i,
                              n.trees = best_iter,
                              i.var = 1,
                              return.grid = TRUE,
                              continuous.resolution = res,
                              type = "response")[,1]
    partDep_2[,1] <- plot.gbm(brt_i,
                              n.trees = best_iter,
                              i.var = 2,
                              return.grid = TRUE,
                              continuous.resolution = res,
                              type="response")[,1]
    partDep_3[,1] <- plot.gbm(brt_i,
                              n.trees = best_iter,
                              i.var = 3,
                              return.grid = TRUE,
                              continuous.resolution = res,
                              type = "response")[,1]
    partDep_4[,1] <- plot.gbm(brt_i,
                              n.trees = best_iter,
                              i.var = 4,
                              return.grid = TRUE,
                              type = "response")[,1]
    TempWarm_AnnDEF_int[,1]  <- plot.gbm(brt_i,
                                         i.var = c(1,2),
                                         n.trees = best_iter,
                                         continuous.resolution = res,
                                         return.grid = TRUE, 
                                         type = "response")[,1]
    TempWarm_AnnDEF_int[,2]  <- plot.gbm(brt_i,
                                         i.var = c(1,2),
                                         n.trees = best_iter,
                                         continuous.resolution = res,
                                         return.grid = TRUE, 
                                         type = "response")[,2]
  }
  
  # FOR ALL ITERATIONS RECORD THE PREDICTED MARGINAL PROBABILITY OF FIRE OCCURRENCE
  # FOR EACH COVARIATE AND THE INTERACTION BETWEEN THE TEMPWARM AND ANNDEF IN THE
  # ALLOCATED COLUMNS.
  partDep_1[,(i+1)] <- plot.gbm(brt_i,
                                n.trees = best_iter,
                                i.var = 1,
                                return.grid = TRUE,
                                continuous.resolution = res,
                                type = "response")[,2]
  partDep_2[,(i+1)] <- plot.gbm(brt_i,
                                n.trees = best_iter,
                                i.var = 2,
                                return.grid = TRUE,
                                continuous.resolution = res,
                                type = "response")[,2]
  partDep_3[,(i+1)] <- plot.gbm(brt_i,
                                n.trees = best_iter,
                                i.var = 3,
                                return.grid = TRUE,
                                type = "response",
                                continuous.resolution = res)[,2]
  partDep_4[,(i+1)] <- plot.gbm(brt_i,
                                n.trees = best_iter,
                                i.var = 4,
                                return.grid = TRUE,
                                type = "response")[,2]
  TempWarm_AnnDEF_int[,i+2]  <- plot.gbm(brt_i,i.var = c(1,2),
                                         n.trees = best_iter,
                                         continuous.resolution = res,
                                         return.grid = TRUE, 
                                         type = "response")[,3]
  
  # OUTPUT PREDICTED PROBABILITY RASTER MAP AS A GEOTIFF FILE
  writeRaster(pred_map,
              paste(save_dir,"pred_map_",i,".tif",sep = ""),
              format = "GTiff",
              overwrite = TRUE)
  
  # OUTPUT BRT AS A RDATA FILE
  save(brt_i,file = paste(save_dir,"brt_",i,".RData",sep = ""))
  
  # PRINT PROGRESS OF ANALYSIS TO COMMAND PROMPT
  cat(paste("BRT ",i," of ",nmodels," ...", sep = ""),"\n")
  if (i == nmodels){
    cat("Finished!")
  }
}

#-----------------------------------------------------------------------------------#
################### EXPORT DATASETS AS CSV FILES ####################################
#-----------------------------------------------------------------------------------#

setwd(save_dir) # CHANGE DIRECTORY TO WHERE THE RESULTS SHOULD
                # BE STORED

# AUC VALUES
colnames(AUC_i) <- c("AUC")
write.csv(x = AUC_i,file = "AUC.csv") 

# RELATIVE INFLUENCE VALUES\
colnames(relInf) <- c(as.character(varinf$var[I]))
write.csv(x = relInf,file = "relInf.csv") 

# NUMBER OF OPTIMUM REGRESSION TREES IN EACH BRT\
colnames(nTrees) <- c("n.trees")
write.csv(x = nTrees,file = "nTrees.csv") 

# CLASSIFICATION ACCURACY RATES
colnames(class_rates) <- c("TPR","FPR","FNR","TNR")
write.csv(x = class_rates,file = "class_rates.csv")

# THRESHOLDS FOR AUC VALUES
colnames(thresh) <- c("probability_threshold")
write.csv(x = thresh,file = "thresholds.csv") 

# TRAINING ERROR INFO FOR EACH BRT
write.csv(x = model_train_err,file = "model_train_err.csv",
                              row.names = FALSE)

# CV VALIDATION ERROR INFO FOR BRTS
write.csv(x = model_valid_err,file = "model_valid_err.csv",
                              row.names = FALSE) 

# OBSERVED AND PREDICTED FRP DATA IN TWO-COLUMNS
colnames(frp_t) <- c("obs","pred")
write.csv(x = frp_t,file = "frp_t.csv")

# OBSERVED FRP DATA ORGANIZED BY ECOREGION
rownames(frp_o) <- ecorinfo[ecor_idx,1]
write.csv(x = frp_o,file = "frp_o.csv") 

# PREDICTED FRP DATA ORGANIZED BY ECOREGION
rownames(frp_p) <- ecorinfo[ecor_idx,1]
write.csv(x = frp_p,file = "frp_p.csv")

# EXPORT PARTIAL DEPENDENCE PLOTS
# PARTIAL DEPENDENCE FOR TEMPWARM
write.csv(x = partDep_1,file = "partDep_TempWarm.csv")
# PARTIAL DEPENDENCE FOR ANNUAL MOISTURE DEFICIT
write.csv(x = partDep_2,file = "partDep_AnnDEF.csv")
# PARTIAL DEPENDENCE FOR TOPOGRAPHIC RUGGEDDNESS
write.csv(x = partDep_3,file = "partDep_TR.csv")
# VEGETATION PARTIAL DEPENDENCE DATA
write.csv(x = partDep_4,file = "partDep_Veg.csv") 
# INTERACTION PARTIAL DEPENDENCE
write.csv(x = TempWarm_AnnDEF_int ,
          file = "TempWarm_AnnDEF_int.csv") 

# RESET WORKING DIRECTORY TO ORIGINAL DIRECTORY
setwd(orgdir)