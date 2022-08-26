% FIG_6.m
% Matlab version: R2012a
%
% This script creates Figure 6 in Young et al. 2016. 
% 
% Fig. 6: Projected fire rotation periods for three different time periods 
% in the 21st century from the AK model. The left-most column represents 
% historical observed (first row) and predicted (second row) fire rotation 
% periods in Alaska, as a reference.  
%
%    Available in '\GBM_Modelling_Analysis_Scripts\AncillaryData\' folder 
%    in 'Figures' nested dataset.
%    ----------------------------------------------------------------------
%    - 'Lat.tif' - latitude information per pixel for plotting maps
%    - 'Lon.tif' - longitude information per pixel for plotting maps
%    - 'obsfrp.tif' - observed 1950-2009 fire rotation periods
%    - 'ecor_info.xlsx' - ecoregion information
%    - 'masks.mat' - spatial masks for each spatial domain
%
%    Available in '\GBM_Modelling_Analysis_Scripts\AncillaryData\Shapefiles\' 
%    folder in 'Figures' nested dataset.
%    ----------------------------------------------------------------------
%    - 'ecor_outlines_for_fig1.shp' - shape file of ecoregions
%    - 'ak_shape_noislands_noSE.shp' - shapefile for mainland Alaska
%
%    Available in '\AK_FINAL_RESULTS\' folder in 'Generalized Boosted 
%    Models and Analysis Scripts' nested dataset.
%    ----------------------------------------------------------------------
%    - 'med_pred_prob.mat' - median predicted probability of fire
%      occurrence
%
%    Available in '\GBM_Modelling_Analysis_Scripts\AncillaryData\Supp_future_projection_files\'
%    in 'Generalized Boosted Models and Analysis Scripts' nested dataset.
%    ----------------------------------------------------------------------
%    - 'med_fut_pred.mat' - median projected probability of fire
%      occurrence for each GCM and future time period. 
%
% DEPENDENCIES:
%    * akaxes.m - function that plots background map of Alaska
%      Available in 'Files4Figs\Figure_Scripts\Functions\' folder in 
%      'Figures' nested dataset.
%    * Functions written by John Iversen and available
%        from the Mathworkds File Exchange:
%        - cbfreeze.m
%        - freezeColors.m 
%        'http://www.mathworks.com/matlabcentral/fileexchange/'
%    * From mapping toolbox:
%        - geotiffread.m
%
% CITATION:
% Young, A. M. et al. Climatic thresholds shape northern high-latitude fire 
% regimes and imply vulnerability to future climate change. 2016.
% Ecography. DOI: 10.1111/ecog.02205
%
% Created by: Adam Young
% Created on: May 2015
% Edited for publication: January 2016
%
% Contact info: Philip E. Higuera, PhD, philip.higuera[at]umontana.edu
%%
% INITIALIZE WORKSPACE
close all;
clear all;
clc;

wdir = '\...\'; % Parent (working) directory where all subfolders from 
                % nested datasets are stored. This may need to be further 
                % modified below depending on how directories are organized 
                % by user.
                
save_dir = [wdir,'Figures\'];

fut_dir = [wdir,'AncillaryData\Supp_future_projection_files\'];

% Addpath to directory to include cbfreeze, freezeColors, and akaxes
% functions to help make current figure
addpath([wdir,'\Scripts\Functions']);

% Load data
cd([wdir,'AncillaryData\']);
load('masks.mat'); mask = masks.akmask;
Lat = geotiffread('Lat.tif');
Lon = geotiffread('Lon.tif');
obs = geotiffread('obsfrp.tif'); obs(obs == -1) = Inf;

cd([wdir,'AncillaryData\Shapefiles']);
ecoshp = shaperead('ecor_outlines_for_fig1.shp', ...
    'UseGeoCoords',true);
shapefile = [pwd,'\ak_shape_noislands_noSE.shp'];

% Load historical predicted probabality of fire occurrence
cd([wdir,'AK_FINAL_RESULTS\']);
load('med_pred_prob.mat');

% Load projected probabality of fire occurrence per pixel
cd(fut_dir);
load('med_fut_pred.mat');

% Dimensions of panels in figure (i.e., 3 by 4)
dims = [3 4];

% Parameters to create figure 
chg = ((6.5*2.54)/min(dims));
dxdy   = 4;
xstart = 0;
ystart = 0;
btwn_xchg = 0.225;
btwn_ychg = 0.6750;

% Latitude and longitude limits for Alaska maps
latPlot = [58 90];
lonPlot = [-176 -140];

myscale = log([50 3200]); % Scale of colormap
% Colormap for Alaskan maps of fire rotation period
cmap = [0.00 0.00 0.00;
        0.50 0.00 0.00;
        1.00 0.00 0.00;
        1.00 0.50 0.00;
        1.00 1.00 0.00;
        1.00 1.00 1.00];

% Information to draw colorbar in Fig. 6
dffxticks = (myscale(2)-myscale(1))/size(cmap,1);
xticks = myscale(1):dffxticks:(myscale(2));
ticks_wanted = {'50', '100', '200', '400', '800', '>1600', ''};

% Labels for Fig. 6
titl = {'Coolest GCM','Warmest GCM','Median 5 GCMs'};
yrs = {'Obs. (1950-2009) ', ...
       '    2010-2039    ', ...
       '    2040-2069    ', ...
       '    2070-2099    ', ...
       'Pred. (1950-2009)',};

% Initialize figure
figure(6); set(gcf, ...
    'Units','Centimeters', ...
    'Position',[28 5 16.9 13.85], ...
    'Color','w');

% Convert thirty year probability of fire occurrence to fire rotation
% period (FRP)
histfrp = 30./med_pred_prob;

% 3-d array to store historical and future FRPs to be plotted. I put the
% maps in this array as it is easier to 
mtx = NaN(size(mask,1),size(mask,2),12);

mtx(:,:,1) = obs .* mask; % Observed 1950-2009 frp
mtx(:,:,5) = histfrp .* mask; % Historical predicted frp

% Calculate median predicted frp of all five GCMs for each 
% future time period 
mind = squeeze(median(med_fut_pred,3));

% MEDIAN 5 GCMS
mtx(:,:,2)  = ((30./(mind(:,:,1))) .* mask); % 2010-2039
mtx(:,:,3)  = ((30./(mind(:,:,2))) .* mask); % 2040-2069
mtx(:,:,4)  = ((30./(mind(:,:,3))) .* mask); % 2070-2099

% WARMEST GCM
mtx(:,:,6)  = ((30./(med_fut_pred(:,:,2,1))) .* mask); % 2010-2039
mtx(:,:,7)  = ((30./(med_fut_pred(:,:,2,2))) .* mask); % 2040-2069
mtx(:,:,8)  = ((30./(med_fut_pred(:,:,2,3))) .* mask); % 2070-2099

% COOLEST GCM
mtx(:,:,10) = ((30./(med_fut_pred(:,:,5,1))) .* mask); % 2010-2039
mtx(:,:,11) = ((30./(med_fut_pred(:,:,5,2))) .* mask); % 2040-2069
mtx(:,:,12) = ((30./(med_fut_pred(:,:,5,3))) .* mask); % 2070-2099
        
% TRANSFORM PREDICTED FRP/MFI VALUES ON A NATURAL LOGARITHM SCALE
mtx = log(mtx);

% Index for plotting
k = 1;

for r = (dims(1)-1):-1:0 % for each row in figure
    for c = 0:(dims(2)-1) % for each column in figure
        
        % Do not plot anything in lower left hand corner
        if r == 0 && c == 0
            k = k + 1;
            continue;
        end
        
        if c > 0
            xtra_chg = 0.225;
        else
            xtra_chg = 0.0;
        end
        
        % Set up axes for Figure 6
        [~,a1] = akaxes([xstart+(c*dxdy)+(c*btwn_xchg)+xtra_chg ystart+(r*dxdy)+(r*btwn_ychg) dxdy], ...
                        latPlot,lonPlot, ...
                        shapefile);
        % Format axes            
        set(gca, ...
            'YColor','w', ...
            'XColor','w', ...
            'Units','centimeters', ...
            'Color','none');
        
        % Plot FRP map
        h1 = surfm(Lat,Lon,mtx(:,:,k));
        % Set up colormap
        colormap(cmap);
        caxis(myscale);
        
        % Plot ecoregion shapefiles
        g1 = geoshow([ecoshp.Lat],[ecoshp.Lon], ...
            'Color','k', ...
            'LineWidth',0.5);
        hold on;
        
        if c == 1
            text(-0.09,1.15,char(titl(r+1)), ...
                'FontName','Arial', ...
                'FontSize',10, ...
                'Rotation',90, ...
                'HorizontalAlignment','Center', ...
                'BackGroundColor','none');
        end
        
        % Add text and further format axes
        if r == 2
            
            text(-0.035,1.2575,char(yrs(k)), ...
                'FontName','Arial', ...
                'FontSize',10);
        end
        
        if r == 1 && c == 0
            
            % Add colorbar
            h = colorbar;
            set(h, ...
                'Position',[0.07 0.03 0.025 0.25], ...
                'YTick',xticks, ...
                'YTickLabel',ticks_wanted, ...
                'FontName','Arial', ...
                'YDir','reverse', ...
                'XGrid','on', ...
                'YGrid','on', ...
                'ZGrid','on', ...
                'GridLineStyle','-', ...
                'LineWidth',0.72);
            
            text(-0.028,0.85,'FRP (yrs)', ...
                'FontName','Arial', ...
                'Rotation',90);
            
            cbfreeze; % Freeze colorbar
            
            % Reset current panel position
            set(a1,'Position', ...
                [xstart+(c*dxdy)+(c*btwn_xchg)+xtra_chg ystart+(r*dxdy)+(r*btwn_ychg) dxdy dxdy]);
            
            
            text(-0.035,1.2575,char(yrs(k)), ...
                'FontName','Arial', ...
                'FontSize',10);
        end
        
        % Freeze current colors
        freezeColors;
        
        % Increase mapping index
        k = k + 1;
        
    end
end

% SAVE FIGURE 6
cd(save_dir);
set(gcf,'PaperType','usletter','PaperPositionMode','auto');
print('FIG_6','-dtiff','-r600');
print('FIG_6','-djpeg','-r600');
saveas(gcf,'FIG_6.fig');