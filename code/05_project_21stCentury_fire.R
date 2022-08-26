# project_21stCentury_fire.R
# 
# THIS SCRIPT TAKES THE BRT MODELS CONSTRUCTED FROM HISTORICAL TRAINING DATA AND USES
# THEM TO PROVIDE PROJECTIONS OF 21ST-CENTURY FOR THE PROBABILITY OF FIRE OCCURRENCE 
# USING AR5 GCM CLIMATE DATA UNDER THE RCP 6.0 SCENARIO. THE SCRIPT EXPORTS THE 
# PROJECTIONS AS GeoTiffs. R version 3.2.2 (2015-08-14) WAS USED TO CONDUCT THIS 
# ANALYSIS.
#                                                                                   
# FILE REQUIREMENTS:
#
#   Available in the "Downscaled climate data for Alaska (2010-2099)" nested dataset.
#   ---------------------------------------------------------------------------------
#   (1) PROJECTED 21ST-CENTURY CLIMATOLOGIES FOR THREE DIFFERENT TIME PERIODS: 
#       2010-2039, 2040-2069, 2070-2099 AND TWO DIFFERENT CLIMATE VARIABLES: 
#       TempWarm & AnnDEF. THESE DATA ARE IN A GEOTIFF FORMAT. 
#   ---------------------------------------------------------------------------------
#
#   Available in the "\1_AK_FINAL_RESULTS\" folder in the Generalized Boosted Models 
#   and Analysis Scripts" nested dataset.
#   ---------------------------------------------------------------------------------
#   (2) BRT MODELS CREATED USING THE gbm PACKAGE. SINCE WE ONLY PROVIDE PROJECTIONS 
#       FROM THE AK MODEL, WE ONLY READ IN THOSE BRTS CREATED FOR THE AK SPATIAL DOMAIN.
#       THESE FILES ARE LOCATED IN: \1_AK_FINAL_RESULTS\
#       THESE FILES ARE IN AN RDATA FORMAT.
#   ---------------------------------------------------------------------------------
#
# DEPENDENCIES:
#   * gbm package (version 2.1.1)
#   * raster package (version 2.4-20)
#   * foreach package (version 1.4.3)
#   * doSNOW package (version 1.0.14)
#
# CITATION:
# Young, A. M. et al. Climatic thresholds shape northern high-latitude fire 
#     regimes and imply vulnerability to future climate change. 2017. 
#     Ecography 40:606-617. doi:10.1111/ecog.02205
#
# Author(s): Adam Young                                                             
# Date created: January 2015                                                              
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

save_dir <- paste(wdir,"FUTURE_PROJECTIONS/",sep="")
                  
# CREATE DIRECTORY, IF IT DOES NOT ALREADY EXIST, TO STORE RESULTS
direxist <- dir.exists(save_dir)
if (direxist == FALSE){
  dir.create(save_dir,
             recursive = TRUE)
}

# LOAD REQUIRED LIBRARIES AND DATASETS

require(gbm) # gbm: GENERALIZED BOOSTING MODELS
require(raster) # raster PACKAGE
# THE FOLLOWING PACKAGES ALLOW FOR THE USE OF ALL FOUR CORES IN PROCESSOR. THIS 
# IMPROVES THE SPEED OF THE TASKS 
require(foreach)
require(doSNOW)

# LOAD STATIC (I.E., UNCHANGING) VEGETATION AND TR DATA TO DRIVE BRTS. WILL HAVE TO 
# MODIFY THE DIRECTORY NAME TO MATCH THE LOCATION OF THESE DATA.
setwd(paste(wdir,"VegLandscapeData/",
            sep = ""))
mapInfo <- raster("AK_VEG.tif")
akVeg   <- getValues(mapInfo) # [categorical, 1-5] EXPLANATORY
TR      <- getValues(raster("TR.tif")) # [meters] EXPLANATORY

# TIME PERIODS TO PROVIDE PROJECTIONS FOR
yrs = c("2010_2039","2040_2069","2070_2099")

# FIVE GCMS
gcms = c("CCSM4","GFDL-CM3","GISS-E2-R","IPSL-CM5A-LR","MRI-CGCM3")

for (g in 1:length(gcms)){ # FOR EACH GCM ...
  
  # IF FOLDER TO STORE FUTURE PROJECTIONS OF FIRE PROBABILITY DOES NOT CURRENTLY
  # EXIST, CREATE THE FOLDER USING mkdirs 
  if (file.exists(paste(save_dir,gcms[g],sep = "")) == FALSE){
    mkdirs(paste(save_dir,gcms[g],sep = ""))
  }
  
  for (y in 1:length(yrs)){ # FOR EACH OF THE THREE TIME PERIODS ...
    
    # CHANGE WORKING DIRECTORY
    setwd(paste(wdir,"FutureClimatologies/",gcms[g],
                sep = ""))
    # LOAD PROJECTED MOISTURE BALANCE AND TEMPERATURE OF THE WARMEST MONTH DATA
    AnnDEF   <- getValues(raster(paste("AnnDEF_",yrs[y],".tif",sep = "")))
    TempWarm <- getValues(raster(paste("TempWarm_",yrs[y],".tif",sep = "")))
    
    # CREATE DATA FRAME TO CREATE PROJECTED PROBABILITY MAP
    proj_data  <- data.frame(TempWarm = TempWarm,
                             AnnDEF   = AnnDEF,
                             TR       = TR,
                             veg      = factor(akVeg))
    
    # LOAD THE OPTIMAL NUMBER OF REGRESSION TREES INFORMATION FOR  BRTS
    setwd(paste(wdir,"1_AK_FINAL_RESULTS/",sep = ""))
    nTrees <- read.csv("nTrees.csv")[,2] 
    
    # INITIALIZE NUMBER OF CORES TO BE USED IN PARALLEL PROCESSING
    cl <- makeCluster(4)
    registerDoSNOW(cl)
    
    # USE foreach FUNCTION TO RUN PARALLEL PROCESSING TO CREATE FUTURE PROJECTIONS 
    # OF FIRE PROBABILITY FOR EACH GCM
    foreach(i = 1:100,.packages=c("raster","gbm")) %dopar% {
      
      # LOAD BOOSTED REGRESSION TREE MODEL
      load(paste("brt_",i,".RData",sep = ""))
      nTrees_i  <- nTrees[i] # NUMBER OF OPTIMAL REGRESSION TREES FOR CURRENT BRTS
      # CREATE PROJECTED PROBABILTIY MAP USING predict.gbm FUNCTION
      pred_map  <- predict.gbm(brt_i,proj_data,
                               n.trees = nTrees_i,
                               type = "response")
      # CONVERT PREDICTED PROBABILTIY MAP TO RASTER FILE
      pred_map  <- raster(t(matrix(pred_map,ncol(mapInfo),nrow(mapInfo))),
                          xmn = xmin(mapInfo),
                          xmx = xmax(mapInfo),
                          ymn = ymin(mapInfo),
                          ymx = ymax(mapInfo),
                          crs = mapInfo@crs)
      
      # EXPORT RASTER MAP
      writeRaster(pred_map,
                  paste(save_dir,gcms[g],"/","AK_",gcms[g],
                        "rcp60_pred_map_",yrs[y],"_",i,".tif",sep=""),
                  format="GTiff",
                  overwrite=TRUE)
      
    }
    stopCluster(cl) # STOP USING MULTIPLE CORES FOR PARALLEL COMPUTING
  } 
} 
