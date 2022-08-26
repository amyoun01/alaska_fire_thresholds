% summarize_future_projections.m
% Matlab version: R2012a
%
% This script summarizes the projections of future fire probability made by
% the R script: project_21stCentury_fire.R. Specifically, this script 
% calculates the (1) median predicted probability of fire occurrence for  
% three different time periods in the 21st century and for three GCMS and 
% (2) the ratio between predicted historical and predicted future fire 
% rotation periods per pixel.
% 
% FILE REQUIREMENTS:
%
%    Available in '\AK_FINAL_RESULTS\', '\BOREAL_FINAL_RESULTS\', and 
%    '\TUNDRA_FINAL_RESULTS\' folders in 'Generalized Boosted 
%    Models and Analysis Scripts' nested dataset.
%    ----------------------------------------------------------------------
%    (1) Historical predicted probability of fire maps for Alaska
%    (1950-2009) from the AK model. These files are located in the
%    folder: '\AK_FINAL_RESULTS'
%    ----------------------------------------------------------------------
%
%    Available in the following zip drives in the 'Future fire projections 
%    for 2010-2099, Alaska' nested dataset:
%    * 'AK_CCSM4_rcp60_pred_maps_2010_2099.zip'
%    * 'AK_GFDL-CM3_rcp60_pred_maps_2010_2099.zip'
%    * 'AK_GISS-E2-R_rcp60_pred_maps_2010_2099.zip'
%    * 'AK_IPSL-CM5A-LR_rcp60_pred_maps_2010_2099.zip'
%    * 'AK_MRI-CGCM3_rcp60_pred_maps_2010_2099.zip'
%    ----------------------------------------------------------------------    
%    (2) Future projected probabilities of fire occurrence for each GCM,
%    time period, and BRT model, located in:
%    
%    - E.g., 'AK_CCSM4_rcp60_pred_map_2010_2039_1.tif'
%    ----------------------------------------------------------------------
%
%    Available in the '\AncillaryData\' folder in the 'Generalized Boosted 
%    Models and Analysis Scripts' nested dataset.
%    ----------------------------------------------------------------------
%    (3) 'TempWarm_1950_2009.tif'
%    ----------------------------------------------------------------------
%
% DEPENDENCIES:
%     * Function from Matlab mapping toolbox:
%       - geotiffread.m 
%
% CITATION:
% Young, A. M. et al. Climatic thresholds shape northern high-latitude fire 
% regimes and imply vulnerability to future climate change. 2016.
% Ecography. DOI: 10.1111/ecog.02205
%
% Created by: Adam Young
% Created on: January 2015
% Edited for publication: January 2016
%
% Contact info: Philip E. Higuera, PhD, philip.higuera[at]umontana.edu
%%
% Initialize workspace
clear all; % clear workspace
close all; % close current figures
clc; % clear command prompt

% Initialize directory names
wdir = '\...\'; % Parent (working) directory where all subfolders from 
                % nested datasets are stored. This may need to be further 
                % modified below depending on how directories are organized 
                % by user.

% LOAD AK MASK MAP TO GET DIMENSIONS OF 2-KM RASTER OF ALASKA. 
tifInfo = geotiffinfo([wdir,'AncillaryData\TempWarm_1950_2009.tif']);

% SET MAIN DIRECTORIES FOR COLLECTING DATA
main_hist_dir = [wdir,'AK_FINAL_RESULTS\']; % HISTORICAL
main_fut_dir  = [wdir,'\...\']; % Directory where all the zip drives are 
                                % extracted. Needs to be filled in.

% CELL ARRAYS TO ITERATIVELY CHANGE DIRECTORY AND RETRIEVE NEEDED RASTER
% MAPS OF ALASKA.
yrs  = {'2010_2039','2040_2069','2070_2099'};
gcms = {'CCSM4','GFDL-CM3','GISS-E2-R','IPSL-CM5A-LR','MRI-CGCM3'};

% FUTURE PREDICTED: THIS IS A 725 x 687 x 5 x 3 ARRAY. THE FIRST
% DIMENSION AND SECOND DIMENSION ARE FOR THE ROWS AND COLUMNS OF THE RASTER
% OF ALASKA, RESPECTIVELY. THE THIRD DIMENSION REPRESENTS A GIVEN GCM, THE
% FOURTH DIMENSION IS REPRESENTATIVE OF THE FUTURE TIME PERIOD
med_fut_pred = NaN(tifInfo.Height,tifInfo.Width, ... 
                   length(gcms), ...
                   length(yrs));

% RATIOS COMPARING PREDICTED FUTURE AND HISTORIC PROBABILITIES: THIS 
% IS THE RATIO OF FRP_FUTURE/FRP_HISTORIC. SINCE THE FRP
% IS CALCULATED BY TAKING THE INVERSE OF THE ANNUAL PREDICTED PROBABILITY OF
% FIRE OCCURRENCE (I.E. 30/PROB), WE SIMPLY CALCULATE:
%
%                           HISTORIC_PROB/FUTURE_PROB
% 
% AS THIS IS EQUAL TO:
%
%                          (30/FUTURE_PROB)/(30/HISTORIC_PROB)
%                        = (30/FUTURE_PROB)*(HISTORIC_PROB/30)
%                        =  HISTORIC_PROB/FUTURE_PROB

% EMPTY ARRAY TO STORE RATIO RESULTS
ratio = NaN(tifInfo.Height,tifInfo.Width, ...
            length(gcms), ...
            length(yrs));

% ARRAY TO STORE HISTORICAL PREDICTIONS OF THE PROBABILITY OF FIRE
% OCCURRENCE FOR ALL ONE-HUNDRED BRTS
hist_pred  = NaN(tifInfo.Height,tifInfo.Width,100);

% ARRAY TO STORE FUTURE PREDICTIONS OF THE PROBABILITY OF FIRE
% OCCURRENCE FOR ALL ONE-HUNDRED BRTS, FIVE GCMS, AND THREE TIME PERIODS
fut_pred = NaN(tifInfo.Height,tifInfo.Width, ...
    100, ...
    length(gcms), ...
    length(yrs));

% ARRAY TO STORE FUTURE RATIOS OF THE RELATIVE CHANGE IN THE FIRE ROTATION 
% PERIOD FOR ALL ONE-HUNDRED BRTS, FIVE GCMS, AND THREE TIME PERIODS
ratio_i  = NaN(tifInfo.Height,tifInfo.Width, ...
    100, ...
    length(gcms), ...
    length(yrs));

for i = 1:100 % FOR EACH BRT
    cd(main_hist_dir); % CHANGE DIRECTORY 
    hist_files = dir('*tif');
    % LOAD HISTORICAL PROBABILITY MAP
    hist_pred(:,:,i)  = geotiffread(['pred_map_',num2str(i),'.tif']);
    
    for g = 1:length(gcms) % FOR EACH GCM
        
        cd([main_fut_dir,char(gcms(g))]); % CHANGE TO THE DIRECTORY WHERE
                                          % FUTURE PROJECTIONS ARE LOCATED.
                                          % Folder names will need to be
                                          % modified to only label the GCM
                                          % name to use this script in its
                                          % current form.

        for y = 1:length(yrs) % FOR EACH TIME PERIOD
            
            % SET NAME OF CURRENT FILE AND AND LOAD PROJECTED PROBABILTIY 
            %  MAP
            filename = sprintf('AK_%s_rcp60_pred_map_%s_%d.tif', ...
                               char(gcms(g)),char(yrs(y)),i);
            fut_pred(:,:,i,g,y) = geotiffread(filename); % STORE PROJECTED
                                                         % PROBABILITY MAP
            % CALCULATE AND STORE RATIO
            ratio_i(:,:,i,g,y) =   hist_pred(:,:,i) ./ fut_pred(:,:,i,g,y);
        end
        
        
    end
end

% SUMMARIZE PROJECTIONS AND RATIO VALUES USING THE MEDIAN PER PIXEL
for g = 1:length(gcms) % FOR EACH GCM
    for y = 1:length(yrs) % FOR EACH TIME PERIOD     
        med_fut_pred(:,:,g,y) = squeeze(median(fut_pred(:,:,:,g,y),3));
        ratio(:,:,g,y) = squeeze(median(ratio_i(:,:,:,g,y),3));
    end
end
% SAVE OUTPUT AS .mat FILES
cd(['\...\']); % Fill in directory to save these files to.
save('med_fut_pred.mat','med_fut_pred');
save('ratio.mat','ratio');
