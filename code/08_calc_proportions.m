% calc_proportions.m
% Matlab version: R2012a
%
% This script calculates the areal percentages expected to experience a
% given level of relative change in the fire rotation period (FRP) relative
% to the historical (1950-2009) period. This is done for all the median
% projected probability of all five GCMs per pixel, as well as the 
% projections for the coolest and warmest GCMs, for three time periods in 
% the 21st-century (2010-2039,2040-2069,2070-2099). 
% 
% FILE REQUIREMENTS:
%
%    Available in '\AncillaryData\' folder in 'Generalized Boosted Models 
%    and Analysis Scripts' nested dataset.
%    ----------------------------------------------------------------------
%    (1) Spatial masks to remove pixels not included in the ananalysis
%        - 'masks.mat'
%    ----------------------------------------------------------------------
%
%    Available in '\AncillaryData\Supp_future_projection_files\'
%    in 'Generalized Boosted Models and Analysis Scripts' nested dataset.
%    ----------------------------------------------------------------------
%    (2) FRP_historical / FRP_future ratios for each pixel and time period
%        - 'ratio.mat'
%    ----------------------------------------------------------------------
%
% CITATION:
% Young, A. M. et al. Climatic thresholds shape northern high-latitude fire 
% regimes and imply vulnerability to future climate change. 2016.
% Ecography. DOI: 10.1111/ecog.02205
%
% Created by: Adam Young
% Created on: February 2015
% Edited for publication: January 2016
%
% Contact info: Philip E. Higuera, PhD, philip.higuera[at]umontana.edu
%%
% Initialize workspace
clear all; % clear workspace
close all; % close current figures
clc; % clear command prompt

wdir = '\...\'; % Parent (working) directory where all subfolders from 
                % nested datasets are stored. Labeled as [WDIR] in 
                % remainder of script. This may need to be further modified
                % below depending on how directories are organized by user.

% Directory of future projections from BRT
fut_dir = [wdir,'Files4Figs\Supp_future_projection_files\'];
save_dir = '\...\'; % Directory to save results to. Needs to be specified. 

% INITIALIZE VARIABLES
ngcms = 3; % Number of "gcms" (median of all five, coolest, warmest)
ntmp  = 3; % Number of time periods
% Discrete categories for proportions of relative change in FRPs
props = [0.25:0.25:1.75];

% Load data
cd([wdir,'AncillaryData\']);
load('masks.mat');
mask = masks.akmask; N = nansum(mask(:));
cd(fut_dir);
load('ratio.mat');

% Calculate the median predicted probability of fire occurence per pixel
% for each time period in the 21st century.
medrat = squeeze(median(ratio,3));

% Allocate space to store percentage data
vals = NaN*ones(length(props)+1,ntmp*ngcms);

k = 1; % Counting variable to help fill matrix of proportion values
for g = 1:ngcms % for each of 3 gcms
    for t = 1:ntmp % for each of 3 time periods
        if g == 1 % if g == 1, then use the median predicted probability
            mtx = medrat(:,:,t).*mask;
        elseif g == 2 % if g == 2, then use the predicted probability from
                      % the warmest gcm
            mtx = ratio(:,:,2,t).*mask;
        elseif g == 3 % if g == 2, then use the predicted probability from
                      % the coolest gcm
            mtx = ratio(:,:,5,t).*mask;
        end
        for j = 1:(length(props)+1) % for each proportional category
            if j == 1
                % This code finds how many pixels in the current map are in
                % a given proportional (or percentage) class
                vals(j,k) = length(mtx(mtx <= props(j)))/N;
                disp([props(j)])
            elseif j > 1 && j < (length(props)+1)
                vals(j,k) = length(mtx(mtx > props(j-1) & mtx <= props(j)))/N;
                disp([props(j-1) props(j)])
            elseif j == length(props)+1
                vals(j,k) = length(mtx(mtx > props(j-1)))/N;
                disp([props(j-1)])
            end 
        end
        k = k + 1; 
    end
end

% Create column and row labels for exporting data table to excel.
collabels1 = {'Med. 5 GCMS','','', ...
              'Warmest GCM','','', ...
              'Coolest GCM','',''};
collabels2 = {'2010_2039','2040_2069','2070_2099', ...
              '2010_2039','2040_2069','2070_2099', ...
              '2010_2039','2040_2069','2070_2099'};
rowlabels  = {'<=0.25','>0.25 & <=0.50','>0.50 & <=0.75','>0.75 & <=1.00', ...
             '>1.00 & <=1.25','>1.25 & <=1.50','>1.50 & <=1.75','>1.75'};

cd(save_dir);
% Export data as an excel file
xlswrite('prop_vals.xlsx',vals,'Sheet1','B3');
xlswrite('prop_vals.xlsx',collabels1,'Sheet1','B1');
xlswrite('prop_vals.xlsx',collabels2,'Sheet1','B2');
xlswrite('prop_vals.xlsx',rowlabels','Sheet1','A3:A10');
