% corr_results.m
% Matlab version: R2012a
%
% This code serves two purposes: 
% (1) Summarize the median predicted probability of fire occurrence from
% the BRT analysis for the AK, BOREAL, and TUNDRA domains for each 2-km
% pixel.
%
% (2) Calculates Pearson correlations between predicted and observed FRP in 
% AK, and export the observed and predicted values, and these correlations 
% in a .mat file to make figures. 
% 
% FILE REQUIREMENTS:
%    
%    Available in '\AK_FINAL_RESULTS\', '\BOREAL_FINAL_RESULTS\', and 
%    '\TUNDRA_FINAL_RESULTS\' folders in 'Generalized Boosted 
%    Models and Analysis Scripts' nested dataset.
%    ----------------------------------------------------------------------
%    (1) Historical predicted probability of fire maps for Alaska
%    (1950-2009). These files are located in each the following folders for
%    each spatial domain (AK, BOREAL, TUNDRA): 
% 
%     - \AK_FINAL_RESULTS
%     - \BOREAL_FINAL_RESULTS
%     - \TUNDRA_FINAL_RESULTS
%        
%    (2) Predicted and observed FRP results from the BRT analysis. For each
%    output directory of the BRT results we use the following files:
%    
%    'frp_t.csv' - two colum vector of observed and predicted frps
%    'frp_o.csv' - n by 100 matrix of observed FRP for each model (i.e.,
%    AK, BOREAL, and TUNDRA). n = the number of ecoregions in each spatial
%    domain.
%    'frp_p.csv' - n by 100 matrix of predicted FRP for each model (i.e.,
%    AK, BOREAL, and TUNDRA). n = the number of ecoregions in each spatial
%    domain.
%        
% DEPENDENCIES:
%     * Functions from Matlab statistics toolbox:
%         - corr.m
%         - nanmedian.m
%     * Function from Matlab mapping toolboxe:
%         - geotiffread.m 
%
% CITATION:
% Young, A. M. et al. Climatic thresholds shape northern high-latitude fire 
% regimes and imply vulnerability to future climate change. 2016.
% Ecography. DOI: 10.1111/ecog.02205
%
% Created by: Adam Young
% Created on: June 2015
% Edited for publication: January 2016
%
% Contact info: Philip E. Higuera, PhD, philip.higuera[at]umontana.edu
%% INITIALIZE WORKSPACE
clear all; % clear workspace of variables
clc; % clear command prompt
close all; % close all current figure windows

% WORKING DIRECTORY
wdir = '\...\'; % Parent (working) directory where all subfolders from 
                % nested datasets are stored. This may need to be further 
                % modified below depending on how directories are organized 
                % by user.
                
% DIRECTORIES OF THE RESULTS FOR EACH OF THE MODELS (AK, BOREAL, & TUNDRA)
models = {'AK_FINAL_RESULTS', ...
          'BOREAL_FINAL_RESULTS', ...
          'TUNDRA_FINAL_RESULTS'};

nmaps = 100; % NUMBER OF PREDICTED PROBABILITY MAPS

for m = 1:length(models) % FOR EACH MODEL ...
    chgdir = [wdir,char(models(m))];
    % CHANGE TO WORKING DIRECTORY
    cd(chgdir);
    
    % LOAD IN EACH PREDICTED PROBABILITY MAP OF ALASKA. STORE IN A 3-D ARRAY.
    pred_prob = NaN(725,687,nmaps);
    for i = 1:nmaps
        fname = sprintf('pred_map_%d.tif',i);
        pred_prob(:,:,i) = double(geotiffread(fname));
    end
    med_pred_prob = nanmedian(pred_prob,3); % CALCULATE MEDIAN PREDICTED
                                            % PROBABILITY OF FIRE OCCURRING IN
                                            % THIRTY YEARS FOR EACH PIXEL IN
                                            % STUDY AREA FROM 100 PREDICTIONS.
    
    % SAVE NEW MEDIAN PREDICTED PROBABILITY MAP IN CURRENT DIRECTORY
    save('med_pred_prob.mat','med_pred_prob');
    
    % LOAD PREDICTED FIRE ROTATION PERIOD FOR EACH ALASKAN ECOREGION IN STUDY
    % AREA IN TWO DIFFERENT FORMATS:
    
    % (1) ALL PREDICTIONS IN THE FORM OF TWO COLUMNS (COL 1 = OBSERVED, COL 2 =
    % PREDICTED). THIS FORMAT IS NAMED 'frp_t.csv', THE '_t' SIGNIFYING USING
    % THE "TOTAL" NUMBER OF OBSERVATIONS AND PREDICTIONS FROM ALL ECOREGIONS.
    % THESE DATA ARE PLOTTTED IN FIGURE 2. 
    frp_all = xlsread('frp_t.csv',1,'B2:C2201');
    frp_all(frp_all == -9999) = NaN;
    
    % (2) THE SECOND FORMAT IS USED TO BETTER CALCULATE THE MEDIAN OBSERVED AND
    % PREDICTED FRP OF EACH ECOREGION BY FORMATTING THE MATRIX INTO WHERE EACH
    % ROW CONTAINS THE FRP OBSERVATIONS AND PREDICTIONS FOR EACH OF THE 100
    % MODELS, REPRESENTED BY COLUMNS.
    frp_obs  = xlsread('frp_o.csv',1,'B1:CW23'); % OBSERVED FRP
    frp_obs(frp_obs == -9999) = NaN;
    frp_pred = xlsread('frp_p.csv',1,'B1:CW23'); % PREDICTED FRP
    frp_pred(frp_pred == -9999) = NaN;
    
    % ALLOCATION OF SPACE TO STORE PEARSON CORRELATION VALUES
    rhoval = NaN(size(frp_obs,2),1); 
    for i = 1:size(frp_obs,2)
        idx = find(isnan(frp_obs(:,i)) == false); % REMOVE MISSING OBSERVATIONS
        rhoval(i) = corr(frp_obs(idx,i),frp_pred(idx,i));
    end
    
    % SUMMARIZE OBSERVED AND PREDICTED FRP VALUES USING THE MEDIAN
    medobs = nanmedian(frp_obs,2);
    medpred = nanmedian(frp_pred,2);
    
    % STORE RESULTS IN A STRUCTURE TO EXPORT AS A .mat FILE.
    corresults.frp_all = frp_all;
    corresults.medobs = medobs;
    corresults.medpred = medpred;
    corresults.rhoval = rhoval;
    
    save('corresults.mat','corresults'); % SAVE DATA STRUCTURE IN CURRENT
                                         % WORKING DIRECTORY.
end
