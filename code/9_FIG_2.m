% FIG_2.m
% Matlab version: R2012a
%
% This script creates Figure 2 in Young et al. 2016. 
% 
% Fig. 2: Depictions of model performance, including (a) observed fire 
% rotation periods (FRPs) from 1950-2009 for Alaskan ecoregions as a 
% reference, (b-d) model predicted fire rotation periods for each 2 × 2 km 
% pixel in Alaska, and (e-g) plots comparing observed fire rotation periods 
% against model predictions per ecoregion. Grey colored points in panels 
% e-g are individual predictions and observations from the 100 boosted 
% regression tree models (BRTs), while the filled darker colored circles 
% and triangles are the median predicted and observed FRPs from the 100 
% BRTs for boreal and tundra ecoregions, respectively. Pearson correlation 
% coefficients (median) are the median recorded Pearson correlation 
% coefficient from a distribution of 100 Pearson correlations comparing the 
% linear relationship between predicted and observed FRPs for each BRT. 
% The x- and y-axes in panels e-g are on the loge scale. 
% Correlations were calculated on untransformed data. 
%
% FILE REQUIREMENTS:
%
%    Available in '\GBM_Modelling_Analysis_Scripts\AncillaryData\' folder 
%    in 'GBM Modelling and Analysis' nested dataset.
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
%    Available in '\AK_FINAL_RESULTS\', '\BOREAL_FINAL_RESULTS\', and 
%    '\TUNDRA_FINAL_RESULTS\' folders in 'Generalized Boosted 
%    Models and Analysis Scripts' nested dataset.
%    ----------------------------------------------------------------------
%    - 'med_pred_prob.mat' - median predicted probability of fire
%    occurrence for each model - to be plotted
%    - 'corresults.mat' - results from Pearson correlations calculated
%    using the script 'corr_results.m'
%
% DEPENDENCIES:
%       * akaxes.m - function that plots background map of Alaska. 
%         Available in 'Files4Figs\Figure_Scripts\Functions\' folder in 
%         'Figures' nested dataset.
%       * freezeColors.m - function written by John Iversen and available
%         from the Mathworkds File Exchange:
%         'http://www.mathworks.com/matlabcentral/fileexchange/7943-freezecolors---unfreezecolors'
%       * From mapping toolbox:
%          - geotiffread.m
%       * From statistics toolbox:
%          - nanmedian.m
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
% Initialize workspace
clear all; % clear workspace
close all; % close current figures
clc;       % clear command prompt

wdir = '\...\'; % Parent (working) directory where all subfolders from 
                % nested datasets are stored. This may need to be further 
                % modified below depending on how directories are organized 
                % by user.

% LOAD DATA
cd([wdir,'GBM_Modelling_Analysis_Scripts\AncillaryData\']);
Lat     = geotiffread('Lat.tif'); % Latitude data
Lon     = geotiffread('Lon.tif'); % Longitude data
obsfrp = geotiffread('obsfrp.tif');  % observed 1950-2009 FRP
         obsfrp(obsfrp == -9999) = NaN;
         obsfrp(obsfrp == -1) = Inf;

% load ecoregion information
[ecorinfo,~,~] = xlsread('ecor_info.xlsx');

% load spatial masks
load('masks.mat');
akmask = masks.akmask;
boreal = masks.boreal;
tundra = masks.tundra;

% load and name shapefiles
cd([wdir,'GBM_Modelling_Analysis_Scripts\AncillaryData\Shapefiles']);
ecoshp = shaperead('ecor_outlines_for_fig1.shp', ...
    'UseGeoCoords',true);
shapefile = [pwd,'\ak_shape_noislands_noSE.shp'];

% Load median predicted probabilities per pixel for each model
cd([wdir,'AK_FINAL_RESULTS\']); 
ak_hist = importdata('med_pred_prob.mat');
cd([wdir,'BOREAL_FINAL_RESULTS\']); 
bor_hist = importdata('med_pred_prob.mat');
cd([wdir,'TUNDRA_FINAL_RESULTS\']); 
tun_hist = importdata('med_pred_prob.mat');

% ADD PATH FOR FOLDER THAT CONTAINS akaxes.m FUNCTION
addpath([wdir,'Figure_Scripts\Functions']);

% Directory to export figures to 
save_dir = [wdir,'Figures\'];
direxist = exist(save_dir,'dir');
if direxist == 0
    mkdir(save_dir); % create directory if it doesn't exist
end

% PARAMETERS FOR INITIAL FIGURE SET UP
fig_lim_r = 28;
fig_lim_l = 2;
fig_ychg  = 10.3;
fig_xchg  = 16.9;

% INITIALIZE FIGURE 2
figure(2);
set(gcf,'Units','Centimeters', ...
    'Position',[fig_lim_r fig_lim_l fig_xchg fig_ychg], ...
    'Color',[1.0 1.0 1.0]);

% PARAMETERS TO CONTROL PANELS WITHIN FIGURE 2
dxdy = 4.0; % FOR SQUARE PANELS THE WIDTH AND HEIGHT IN CM
xstart = 0; % X-POSITION OF FIRST PANEL
ystart = 1.2; % Y-POSITION OF FIRST PANEL
btwn_figs_y = 0.5; % HEIGHT SPACE BETWEEN PANELS
btwn_figs_x = 0.225; % WIDTH SPACE BETWEEN PANELS
cbar_width = 0.5; % COLORBAR WIDTH
cbar_height = 4.8; % COLORBAR HEIGHT

cbar_shift = (fig_xchg - (xstart+2*dxdy+btwn_figs_x) - cbar_width)/(4*dxdy);

% LATITUDE AND LONGITUDE LIMITS FOR akaxes.m
latPlot = [58 90];
lonPlot = [-168 -140];

% SCALE AND COLORMAP FOR FRP MAP
myscale = log([50 3200]);
cmap = [0.00 0.00 0.00;
        0.50 0.00 0.00;
        1.00 0.00 0.00;
        1.00 0.50 0.00;
        1.00 1.00 0.00;
        1.00 1.00 1.00];

% TICK INFORMATION FOR COLORBAR INFORMATION
dffxticks = (myscale(2)-myscale(1))/size(cmap,1);
xticks = myscale(1):dffxticks:(myscale(2));
ticks_wanted = {'50', '100', '200', '400', '800', '>1600', ''};

% x- and y-axis limits for plotting purposes
ylims = [50 25000];
xlims = [50 25000];

% FOLDERS WHERE MODEL PERFORMANCE RESULTS ARE STORED
indfolders = {'AK_FINAL_RESULTS', ...
              'BOREAL_FINAL_RESULTS', ...
              'TUNDRA_FINAL_RESULTS'};

for i = 1:length(indfolders) % FOR EACH MODEL
    
    % CHANGE DIRECTORY
    directory = [wdir,char(indfolders(i))];
    cd(directory);
    
    % LOAD CORRELATION RESULTS
    corresults   = load('corresults.mat');
    corresults   = corresults.corresults;
    rhoval = nanmedian(corresults.rhoval);
    
    % USE ECOREGION TYPE INFORMATION TO PLOT DIFFERENT SYMBOLS FOR BOREAL
    % FOREST AND TUNDRA ECOREGIONS
    [~,order] = sort(ecorinfo(:,1));
    type = ecorinfo(order,2);
%     type(type==3) = [];
    
    % LOAD PREDICTED AND OBSERVED FIRE ROTATION PERIOD INFORMATION
    results  = [corresults.medobs,corresults.medpred];
    x = results(:,1);
    y = results(:,2);
    
    data = corresults.frp_all;
    
    % PANEL FOR CURRENT PLOT
    subplot(2,4,i+5);
    
    % FOR EACH  MODEL PLOT PREDICTED VS. OBSERVED FRP DATA
    if i == 1 
        
        bor_x = (x(type==2 | type==3));
        bor_y = (y(type==2 | type==3));
        tun_x = (x(type==1));
        tun_y = (y(type==1));
        
        plot(data(:,2),data(:,1),'.', ...
            'MarkerEdgeColor',[0.7 0.7 0.7]);
        line(xlims, ...
            ylims, ...
            'LineStyle','--', ...
            'Color','k', ...
            'LineWidth',1.44)
        hold on;
        plot(bor_y,bor_x,'ko','LineWidth',1,'MarkerFaceColor','k');
        plot(tun_y,tun_x,'k^','LineWidth',1,'MarkerFaceColor','k');
        
        set(gca,'Units','Centimeters', ...
            'Position',[xstart+dxdy+btwn_figs_x ystart dxdy dxdy]);
        
        text(1.4*min(xlims),0.6*max(ylims),'(e)', ...
            'FontSize',11, ...
            'FontName','Arial',...
            'FontWeight','Bold');
        text(0.0334*max(xlims),0.00334*max(ylims), ...
            sprintf('\\it\\rho\\rm_m_e_d_i_a_n= %1.2f',rhoval), ...
            'FontName','Arial', ...
            'FontSize',9);
        
        ylabel('Obs. FRP (yrs)', ...
            'FontName','Arial', ...
            'FontSize',10);
        
    elseif i == 2
        
        plot(data(:,2),data(:,1),'.','MarkerEdgeColor',[0.7 0.7 0.7]);
        line(xlims, ...
            ylims, ...
            'LineStyle','--', ...
            'Color','k', ...
            'LineWidth',1.44)        
        hold on;
        plot(y,x,'ko','LineWidth',1,'MarkerFaceColor','k'); hold on;
        
        set(gca,'Units','Centimeters',...
            'position',[xstart+2*dxdy+2*btwn_figs_x ystart dxdy dxdy],...
            'YTickLabel',{''});
        axis square;
        
        text(1.4*min(xlims),0.6*max(ylims),'(f)', ...
            'FontSize',11, ...
            'FontName','Arial',...
            'FontWeight','Bold');
        text(0.0334*max(xlims),0.00334*max(ylims), ...
            sprintf('\\it\\rho\\rm_m_e_d_i_a_n= %1.2f',rhoval), ...
            'FontName','Arial', ...
            'FontSize',9);

        xlabel('Pred. FRP (yrs)', ...
            'FontName','Arial', ...
            'FontSize',10);
        
    elseif  i == 3
        
        plot(data(:,2),data(:,1),'.','MarkerEdgeColor',[0.7 0.7 0.7]);
        line(xlims, ...
            ylims, ...
            'LineStyle','--', ...
            'Color','k', ...
            'LineWidth',1.44)        
        hold on;
        plot(y,x,'k^','LineWidth',1,'MarkerFaceColor','k'); hold on;
        
        set(gca,'Units','Centimeters',...
            'Position',[xstart+(3*dxdy)+(3*btwn_figs_x) ystart dxdy dxdy],...
            'YTickLabel',{''});
        
        text(1.4*min(xlims),0.6*max(ylims),'(g)', ...
            'FontSize',11, ...
            'FontName','Arial',...
            'FontWeight','Bold');
        text(0.0334*max(xlims),0.00334*max(ylims), ...
            sprintf('\\it\\rho\\rm_m_e_d_i_a_n= %1.2f',rhoval), ...
            'FontName','Arial', ...
            'FontSize',9);
        
    end
    
    % SET AXES 
    xlim(xlims);
    ylim(ylims);
    
    set(gca,'XScale','log',...
        'YScale','log', ...
        'XTick',[100 1000 10000], ...
        'YTick',[100 1000 10000], ...
        'XMinorTick','off', ...
        'YMinorTick','off', ...
        'LineWidth',0.72);    
end
%%
% BELOW, PLOT EACH MAP OF HISTORICAL PREDICTED PROBABILITIES AND HISTORICAL
% OBSERVED FRP

% AK MAP MAP OF PREDICTED PROBABILITIES
[~,hx2] = akaxes([xstart+dxdy+btwn_figs_x ystart+dxdy+btwn_figs_y dxdy], ...
                 latPlot,lonPlot, ...
                 shapefile);
box off;
axis tight;
set(gca,'Ycolor','w',...
    'Xcolor','w');

surfm(Lat,Lon,log(30./ak_hist).*akmask);
colormap(cmap);
caxis(myscale);

geoshow([ecoshp.Lat],[ecoshp.Lon],'Color','k','LineWidth',0.5);

title('AK', ...
    'FontSize',10, ...
    'FontName', 'Arial')

freezeColors;
textm(70,-170,'(b)', ...
    'FontName','Arial', ...
    'FontWeight','Bold', ...
    'FontSize',11);

hold on;
%%
% BOREAL MAP OF PREDICTED PROBABILITIES
[~,hx3] = akaxes([xstart+2*dxdy+2*btwn_figs_x ystart+dxdy+btwn_figs_y dxdy], ...
                 latPlot,lonPlot, ...
                 shapefile);

box off;
axis tight

set(gca,'YColor','w',...
    'XColor','w');

surfm(Lat,Lon,log(1./(bor_hist./30)).*boreal);
colormap(cmap);
caxis(myscale)

geoshow([ecoshp.Lat],[ecoshp.Lon],'Color','k','LineWidth',0.5);

title('BOREAL', ...
    'FontSize',10, ...
    'FontName', 'Arial')

freezeColors;
textm(70,-170,'(c)','FontName','Arial','FontWeight','Bold','FontSize',11);

% TUNDRA MAP OF PREDICTED PROBABILITIES
[~,hx4] = akaxes([xstart+(3*dxdy)+(3*btwn_figs_x) ...
                  ystart+dxdy+btwn_figs_y dxdy], ...
                 latPlot,lonPlot, ...
                 shapefile);

box off;
axis tight;

set(gca,'YColor','w',...
    'XColor','w');

h2 = surfm(Lat,Lon,log(1./(tun_hist./30)).*tundra);
colormap(cmap);
caxis(myscale)

geoshow([ecoshp.Lat],[ecoshp.Lon],'Color','k','LineWidth',0.5);

title('TUNDRA', ...
    'FontSize',10, ...
    'FontName', 'Arial')

freezeColors;
textm(70,-170,'(d)','FontName','Arial','FontWeight','Bold','FontSize',11);

titl = get(gca,'Title');
set(titl,'Units','Centimeters');
pos_titl = get(titl,'Position');
df_titl_orig = pos_titl(2) - dxdy;

% OBS FRP 
[~,hx1] = akaxes([xstart ...
                  ystart+dxdy+btwn_figs_y ...
                  dxdy], ...
                 latPlot,lonPlot, ...
                 shapefile);
pos = get(gca,'Position');
h1 = surfm(Lat,Lon,log(obsfrp).*akmask);
colormap(cmap);
caxis(myscale);
h=colorbar;
set(gca, ...
    'YColor','w', ...
    'XColor','w', ...
    'Position',pos);
set(h, ...
    'Units','Centimeters', ...
    'Position',[pos(1)+dxdy/4 ystart cbar_width 4], ...
    'YTick',xticks, ...
    'YTickLabel',ticks_wanted, ...
    'FontName','Arial', ...
    'FontSize',9, ...
    'YDir','reverse', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'ZGrid','on', ...
    'GridLineStyle','-', ...
    'LineWidth',0.72);

title('Observed (1950-2009)', ...
    'FontSize',10, ...
    'FontName', 'Arial')

geoshow([ecoshp.Lat],[ecoshp.Lon],'Color','k','LineWidth',.5);

textm(70,-170,'(a)', ...
    'FontName','Arial', ...
    'FontWeight','Bold', ...
    'FontSize',11);

text(-0.075,0.85, ...
    'FRP (yrs)', ...
    'FontName','Arial', ...
    'Rotation',90);

% SAVE FIGURE AS .FIG AND .TIF FILES AT 600 DPI RESOLUTION
cd(save_dir);
set(gcf,'PaperType','usletter','PaperPositionMode','auto');
print('FIG_2','-dtiff','-r600');
print('FIG_2','-djpeg','-r600');
saveas(gcf,'FIG_2.fig');
