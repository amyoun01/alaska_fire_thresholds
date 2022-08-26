% FIG_7.m
% Matlab version: R2012a
%
% This script creates Figure 6 in Young et al. 2016. 
% 
% Fig. 7: Relative change in the fire rotation period 
% (FRPFUTURE / FRPHISTORICAL) per pixel for three different time periods in 
% the 21st century . An increase in the future probability of fire is 
% depicted on a nonlinear scale (e.g., a ratio of 0.5 is equal to a 100% 
% increase in area burned, while a ratio 2 is equal to a 50% decrease). An 
% increase in the future probability of fire is shown by warmer colors and 
% decreasing fire rotation periods (relative difference < 1.0), while a 
% decrease in future fire activity is shown by cooler colors and increasing 
% fire rotation periods (relative difference > 1.0). Pie charts depict the 
% proportions of all pixels in the study domain projected to experience a 
% given level of relative change.
%
%    Available in '\Files4Figs\' folder in 'Figures' nested dataset.
%    ----------------------------------------------------------------------
%    - 'Lat.tif' - latitude information per pixel for plotting maps
%    - 'Lon.tif' - longitude information per pixel for plotting maps
%    - 'obsfrp.tif' - observed 1950-2009 fire rotation periods
%    - 'ecor_info.xlsx' - ecoregion information
%    - 'masks.mat' - spatial masks for each spatial domain
%
%    Available in '\Files4Figs\Shapefiles\' folder in 'Figures' nested 
%    dataset.
%    ----------------------------------------------------------------------
%    - 'ecor_outlines_for_fig1.shp' - shape file of ecoregions
%    - 'ak_shape_noislands_noSE.shp' - shapefile for mainland Alaska
%
%    Available in '\Files4Figs\Supp_future_projection_files\'
%    in 'Figures' nested dataset.
%    ----------------------------------------------------------------------
%    - 'ratio.mat' - ratio between future and historical predicted FRPs per
%    pixel
%
% DEPENDENCIES:
%    * akaxes.m - function that plots background map of Alaska.
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

% Set folder names
wdir = '\...\'; % Parent (working) directory where all subfolders from 
                % nested datasets are stored. This may need to be further 
                % modified below depending on how directories are organized 
                % by user.
                
save_dir = [wdir,'GBM_Modelling_Analysis_Scripts\'];

% These files are the listed in the 
fut_dir = [wdir,'AncillaryData\Supp_future_projection_files\'];

% Addpath to directory to include cbfreeze, freezeColors, and akaxes
% functions to help make current figure
addpath([wdir,'Scripts\Functions']);

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

cd(fut_dir);
load('ratio.mat');
prop_vals = xlsread('prop_vals.xlsx','Sheet1','B3:J10');

% Discrete classes of the relative change in future FRPs relative to the
% past
prop_idx = [0.25:0.25:1.75];

% Three time periods of interest
yrs = {'2010-2039','2040-2069','2070-2099'};

% Parameters to create figure and figure panels
dims = [3 3];
dxdy = 4.0;
xstart = 2.55;
ystart = 2.3;
btwn_xchg = 0.90;
btwn_ychg = 0.45;

% Colormap information for AK ratio maps
cmap = flipud(jet(8));
myscale = [0 2.0];

% Lat and lon limits for AK maps
latPlot = [58 90];
lonPlot = [-176 -140];

titl = {'Coolest GCM','Warmest GCM','Median 5 GCMs'};

xticks = [0.25:0.25:1.75];
ticks_wanted = {'\leq 0.25','0.50','0.75','1.00','1.25','1.50','> 1.75'};

% Initialize figures
figure(7); set(gcf, ...
    'Units','Centimeters', ...
    'Position',[28 5 16.9 15.8], ...
    'Color','w');

mtx = NaN(size(mask,1),size(mask,2),dims(1)*dims(2));

medrat = squeeze(median(ratio,3));

% MEDIAN 5 GCMS
mtx(:,:,1)  = medrat(:,:,1) .* mask; % 2010-2039
mtx(:,:,2)  = medrat(:,:,2) .* mask; % 2040-2069
mtx(:,:,3)  = medrat(:,:,3) .* mask; % 2070-2099

% WARMEST GCM
mtx(:,:,4)  = ratio(:,:,2,1) .* mask; % 2010-2039
mtx(:,:,5)  = ratio(:,:,2,2) .* mask; % 2040-2069
mtx(:,:,6)  = ratio(:,:,2,3) .* mask; % 2070-2099

% COOLEST GCM
mtx(:,:,7) = ratio(:,:,5,1) .* mask; % 2010-2039
mtx(:,:,8) = ratio(:,:,5,2) .* mask; % 2040-2069
mtx(:,:,9) = ratio(:,:,5,3) .* mask; % 2070-2099

k = 1;
for r = (dims(1)-1):-1:0
    
    for c = 0:(dims(2)-1)
        
        [~,a1] = akaxes([xstart+(c*dxdy)+(c*btwn_xchg) ...
                         ystart+(r*dxdy)+(r*btwn_ychg) dxdy], ...
                         latPlot,lonPlot, ...
                         shapefile);
        set(gca, ...
            'YColor','w', ...
            'XColor','w', ...
            'Units','centimeters', ...
            'Color','none');
        h1 = surfm(Lat,Lon,mtx(:,:,k));
        colormap(cmap);
        caxis(myscale);
        freezeColors;
        hold on;
        g1 = geoshow([ecoshp.Lat],[ecoshp.Lon], ...
            'Color','k', ...
            'LineWidth',0.5);
        
        if r == 2
            text(-0.02,1.2575,char(yrs(k)), ...
                'FontName','Arial', ...
                'FontSize',10);
        end
        
        if c == 0
            text(-0.15,1.1,char(titl(r+1)), ...
                'FontName','Arial', ...
                'FontSize',10, ...
                'Rotation',90, ...
                'HorizontalAlignment','Center', ...
                'BackGroundColor','none');
        end
        
        if r == 0 && c == 1
            
            h = colorbar;
            set(h, ...
                'Units','Centimeters', ...
                'Location','SouthOutside', ...
                'XTick',xticks, ...
                'XTickLabel',{''}, ...
                'FontName','Arial', ...
                'FontSize',8, ...
                'YDir','reverse', ...
                'XGrid','on', ...
                'YGrid','on', ...
                'ZGrid','on', ...
                'GridLineStyle','-', ...
                'LineWidth',0.72);
            
            cbar_title = get(h,'Title');
            set(cbar_title, ...
                'String','Relative change in FRP (Future/Historical)', ...
                'FontSize',10, ...
                'FontName','Arial');
            
            set(gca,'Position', ...
                [xstart+(c*dxdy)+(c*btwn_xchg) ...
                 ystart+(r*dxdy)+(r*btwn_ychg) ...
                 dxdy ...
                 dxdy]);
            
            set(h,'Position',[xstart+dxdy-0.1 ystart-2.75*btwn_ychg 6 0.4]);
            
            text_idx = 1;
            for t = -0.25:(6/8):(6-2*(6/8))
                text(t,-1.275,char(ticks_wanted(text_idx)), ...
                    'Rotation',300, ...
                    'Units','Centimeters');
                text_idx = text_idx + 1;
            end
            
        end
        
        axes('Units','Centimeters', ...
            'Position',[xstart+(c*dxdy)+(c*btwn_xchg)-1.3*btwn_xchg ...
                        ystart+(r*dxdy)+(r*btwn_ychg)+6*btwn_ychg ...
                        1.4 ...
                        1.4]);
        
        prop_i = prop_vals(:,k);
        
        prop_i(prop_i<0.00001) = 0.00001;
        pie(prop_i,repmat({''},1,length(prop_idx)+1));
        colormap(cmap)

        k = k + 1;
    end
end

% SAVE FIGURES
cd(save_dir);
set(gcf,'PaperType','usletter','PaperPositionMode','auto');
print('FIG_7','-dtiff','-r600');
print('FIG_7','-djpeg','-r600');
saveas(gcf,'FIG_7.fig');