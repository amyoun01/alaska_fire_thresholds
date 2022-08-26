Young, A. M. et al. 2017. Climatic thresholds shape northern high-latitude fire regimes and imply 
vulnerability to future climate change. Ecography. https://doi.org/10.1111/ecog.02205 

Ecography DOI: 10.1111/ecog.02205
Dryad Repository DOI: 10.5061/dryad.r217r 
NSF Arctic Data Center DOI: 	

Contact: Philip E. Higuera, philip.higuera[at]umontana.edu

########################################################################################################################
These scripts were originally published with the associated datasets on the Dryad Data Repository and this work is 
licensed under a CC0 1.0 Universal (CC0 1.0) Public Domain Dedication license. Please cite both the published paper
in Ecogrpahy and the Dryad Repository when using or referencing the data and code used to complete this work. 

Paper citation: Young, A.M., Higuera, P.E., Duffy, P.A. and Hu, F.S. (2017), Climatic thresholds shape northern 
                       high-latitude fire regimes and imply vulnerability to future climate change. Ecography, 40: 
                       606-617. https://doi.org/10.1111/ecog.02205

 Dataset citation: Young, Adam M.; Higuera, Philip E.; Duffy, Paul A.; Hu, Feng Sheng (2016), Data from: Climatic 
                       thresholds shape northern high-latitude fire regimes and imply vulnerability to future climate 
                       change, Dryad, Dataset, https://doi.org/10.5061/dryad.r217r
########################################################################################################################

######################################################################################################################### 
There are two subsets of files in this project folder: (1) Fourteen R (.R) and Matlab (.m) files sequentially numbered to 
reproduce the results in Young et al. 2017, and (2) a folder containing custom written functions that support the sixteen 
scripts. If the scripts are run in the same order as they are numbered then the exact results and figures from the main 
text will be reproduced. Please note, there are several functions that need to be downloaded from the MATLAB File Exchange
in order to successfully run the scripts used to make several of the figures. Please see the comments in each script for 
details. 
#########################################################################################################################

-------------------------------------------------------------------------------------------------------------------------
Analysis scripts:

	* 1_AK_BRTS.R – Script that uses climatologies created in ‘1_create_hist_ak_climatologies.m’, vegetation data, and 
          topographic ruggedness data in conjunction with the gbm package (v2.1.1) to create the set of 100 BRTs that comprise 
          the AK model.
	* 2_BOREAL_BRTS.R – Script that uses climatologies created in ‘1_create_hist_ak_climatologies.m’, vegetation data, and 
          topographic ruggedness data in conjunction with the gbm package (v2.1.1) to create the set of 100 BRTs that comprise 
          the BOREAL model.
	* 3_TUNDRA_BRTS.R – Script that uses climatologies created in ‘1_create_hist_ak_climatologies.m’, vegetation data, and 
          topographic ruggedness data in conjunction with the gbm package (v2.1.1) to create the set of 100 BRTs that comprise 
          the TUNDRA model.
	* 4_Identify_Climatic_Thresholds.R – Script that uses the segmented package (v0.5-1.4) to identify climatic thresholds 
          from BRT partial dependence results. Partial dependence results were exported as .csv files in the three BRT scripts 
          (i.e., #3-5). 
	* 5_project_21stCentury_fire.R – Script that uses the BRTs that comprise the AK model (created in ‘3_AK_BRTS.R’) to 
          project the 30-yr probability of fire occurrence per pixel for Alaska in the 21st-century. Future climatologies were 
          created using the ‘2_create_future_ak_climatologies.m’ script. 
	* 6_corr_results.m – Main purpose is to calculate the Pearson correlation between predicted and observed fire rotation 
          period (FRP) estimates for Alaskan ecoregions. FRP estimates were exported in the BRT scripts (i.e., #3-5). 
	* 7_summarize_future_projections.m – Takes BRT-generated future projections of the 30-yr probability of fire occurrence 
          and calculates (1) the median predicted 30-yr probability of fire occurrence of all 100 BRT projections for all five GCMs 
          and for three time periods, and (2) the median ratio between future and historical FRPs for all GCMs and time periods.
	* 8_calc_proportions.m – Calculate the percentage of pixels that occurs in eight different discrete classes of relative 
          change in the FRP (i.e., ratio).
	* 11_FIG_2.m – Creates Figure 2 in Young et al. (2016).
	* 12_FIG_3.m – Creates Figure 3 in Young et al. (2016).
	* 13_FIG_4.m – Creates Figure 4 in Young et al. (2016).
	* 14_FIG_5.m – Creates Figure 5 in Young et al. (2016).
	* 15_FIG_6.m – Creates Figure 6 in Young et al. (2016).
	* 16_FIG_7.m – Creates Figure 7 in Young et al. (2016).

-------------------------------------------------------------------------------------------------------------------------
Functions: 

	* akaxes.m – plot background map of Alaska.
	* aucroc.R – Calculates AUC value and classification rates.This function was inspired and based off of the 
          "evaluate" function in the 'dismo' package (Robert J. Hijmans et al. 2016)
	* createClimatologies.m – Function that takes gridded GeoTiff files at a monthly or annual resolution and calculates climatological averages 
          per pixel for a user-defined time period.

        ----------------------------------------------------------------------------------------------------------------------------
        PLEASE NOTE, THE FOLLOWING FUNCTIONS MUST BE DOWNLOADED FROM THE MATLAB FILE EXCHANGE TO SUCCESSFULLY RUN THE FIG_2, FIG_3, 
	FIG_6, AND FIG_7 SCRIPTS.
        ----------------------------------------------------------------------------------------------------------------------------

	* barwitherr.m – Creates bar graph with error bars  - used in Fig. 3. Created by Martina Callaghan and obtained from Matlab 
          file exchange (http://www.mathworks.com/matlabcentral/fileexchange/30639-barwitherr-errors-varargin-). 
	* cbfreeze.m – Freezes colors on colorbar. Created by Carlos Adrian Vargas Aguilera and obtained from Matlab file exchange 
          (http://www.mathworks.com/matlabcentral/fileexchange/24371-colormap-and-colorbar-utilities--jul-2014-/content/cm_and_cb_utilities/cbfreeze.m). 
	* freezeColors.m – Freezes colors on current plot. Created by John Iversen and obtained from Matlab file exchange 
          (http://www.mathworks.com/matlabcentral/fileexchange/7943-freezecolors---unfreezecolors).

-------------------------------------------------------------------------------------------------------------------------
References:

Robert J. Hijmans, Steven Phillips, John Leathwick and Jane Elith (2016). dismo: Species Distribution Modeling. 
       R package version 1.0-15. http://CRAN.R-project.org/package=dismo
