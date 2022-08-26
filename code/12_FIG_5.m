% FIG_5.m
% Matlab version: R2012a
%
% This script creates Figure 5 in Young et al. 2016.
%
% Fig. 5: Interactions between the mean temperature of the warmest month 
% (TWARM) and annual moisture availability (P-PETANN), and the 30-yr probability 
% of fire occurrence per pixel for the (a) AK, (b) BOREAL, and (c) TUNDRA 
% models. The response surface represents the median predicted probability 
% of fire occurrence from 100 boosted regression tree models for each model 
% type. Darker (lighter) colors in the response surface represent higher 
% (lower) probabilities of fire occurrence. A lowess function (span=0.1) 
% was used to smooth the response surface.
%
% FILE REQUIREMENTS:
%
%    Available in '\GBM_Modelling_Analysis_Scripts\AncillaryData\' folder 
%    in 'Figures' nested dataset.
%    ----------------------------------------------------------------------
%    - 'climlims.mat' - 
%
%    Available in '\AK_FINAL_RESULTS\', '\BOREAL_FINAL_RESULTS\', and 
%    '\TUNDRA_FINAL_RESULTS\' folders in 'Generalized Boosted 
%    Models and Analysis Scripts' nested dataset.
%    ----------------------------------------------------------------------
%    - 'TempWarm_AnnDEF_int.csv' - TempWarm and AnnDEF interaction 
%                                  partial dependence results
%
% DEPENDENCIES:
%    * From curve fitting toolbox:
%        - smooth.m
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
wdir = '\...\'; % Parent (working) directory where all subfolders from 
                % nested datasets are stored. This may need to be further 
                % modified below depending on how directories are organized 
                % by user.

% change directory
cd([wdir,'Files4Figs\']);

% Directory of results
save_dir = [wdir,'Figures\'];
modeltype = {'TUNDRA_FINAL_RESULTS', ...
             'BOREAL_FINAL_RESULTS', ...
             'AK_FINAL_RESULTS'};

% axis limits
xlims = [0 20];
ylims = [-500 750];
zlims = [ 0.0 0.08; 0.0 0.6; 0.0 0.5];

% number of ticks per axis
nticks = [4 5 5];


xscl = [0.86 0.98];
yscl = [0.83 0.98];
zscl = [1.3 1.3];

% INITIALIZE FIGURE
figure(5); 
clf; 
set(gcf, ...
    'color','w', ...
    'units','centimeters',...
    'position',[28.0 5.0 8.0 17.2]);

% PARAMETRS TO CREATE FIGURE AND PANELS WITH
ystart = 0.8;
xstart = 1.25;
dx     = 6.0;
dy     = 4.9;
ychg   = 0.65;

% FOR EACH MODEL TYPE
for u = 1:length(modeltype)
    
    % change directory to model folder (e.g., BOREAL)
    chgdir = [wdir,cell2mat(modeltype(u))];
    cd(chgdir);
    
    % Load interaction partial dependence results
    TempWarm_AnnDEF = importdata('TempWarm_AnnDEF_int.csv');
    TempWarm_AnnDEF = TempWarm_AnnDEF.data;
    
    % extract variables from loaded dataset
    TempWarm = TempWarm_AnnDEF(:,1);
    AnnDEF   = TempWarm_AnnDEF(:,2);
    prob     = TempWarm_AnnDEF(:,3:end);
    med_prob = median(prob,2);
    
    % Number of predicted probability points for each covariate
    resolution = size(TempWarm,1)/100;
    
    x = TempWarm(1:resolution);
    y = AnnDEF(1:resolution:((resolution^2)-(resolution-1)));
                    
    Z = reshape(med_prob,resolution,resolution)';
    
    % Meshgrid of x and y data values for plotting
    [X,Y] = meshgrid(x,y);
    
    % Smooth 3-d plane
    [n m] = size(Z);
    smoothInt1 = 0.10;
    smoothInt2 = 0.10;
    Zsmooth = NaN*ones(size(Z));
    for j = 1:m
        Zsmooth(:,j) = smooth(Z(:,j),smoothInt1,'lowess');
        for i = 1:n
            Zsmooth(i,:) = smooth(Z(i,:),smoothInt2,'lowess');
        end
    end
    for j = 1:m
        Zsmooth(:,j) = smooth(Zsmooth(:,j),smoothInt1,'lowess');
        for i = 1:n
            Zsmooth(i,:) = smooth(Zsmooth(i,:),smoothInt2,'lowess');
        end
    end
    
    % Create colormap
    lmt = 0.9;
    lwr = 0.0;
    fireCmap = NaN*ones(resolution,3);
    fireCmap(:,1) = flipud([lwr:(lmt-lwr)/(resolution-1):lmt]');
    fireCmap(:,2) = flipud([lwr:(lmt-lwr)/(resolution-1):lmt]');
    fireCmap(:,3) = flipud([lwr:(lmt-lwr)/(resolution-1):lmt]');
    
    % Set up axes for figures
    v = u - 1; % use to help shift position for each panel
    axes('Units','Centimeters', ...
         'Position',[xstart ystart+v*ychg+v*dy dx dy]);

     idx_frame = [1:5:100,100];
     
    % Create interaction (i.e., surface) plot 
    sf = surf(x(idx_frame),y(idx_frame), ...
              Zsmooth(idx_frame,idx_frame), ...
              'EdgeColor','k', ...
              'FaceAlpha',0.75); 
    hold on;
    
    % Set colormap for surface plot
    colormap(fireCmap)
    
    % SET X-, Y-, AND Z-AXIS LIMITS
    xlim(xlims);
    ylim(ylims);
    zlim(zlims(u,:));
    
    % Add text to figure panels
    if u == 3
        text(xscl(2)*max(xlims),yscl(2)*max(ylims),zscl(2)*max(zlims(u,:)), ...
            '(a)', ...
            'FontName','Arial', ...
            'FontSize',11, ...
            'FontWeight','Bold');
        text(xscl(1)*max(xlims),yscl(1)*max(ylims),zscl(1)*max(zlims(u,:)), ...
            'AK', ...
            'FontName','Arial', ...
            'FontSize',10);
    elseif u == 2
        text(xscl(2)*max(xlims),yscl(2)*max(ylims),zscl(2)*max(zlims(u,:)), ...
            '(b)', ...
            'FontName','Arial', ...
            'FontSize',11, ...
            'FontWeight','Bold');
        text(xscl(1)*max(xlims),yscl(1)*max(ylims),zscl(1)*max(zlims(u,:)), ...
            'BOREAL', ...
            'FontName','Arial', ...
            'FontSize',10);
    elseif u == 1
        text(xscl(2)*max(xlims),yscl(2)*max(ylims),zscl(2)*max(zlims(u,:)), ...
            '(c)', ...
            'FontName','Arial', ...
            'FontSize',11, ...
            'FontWeight','Bold');
        text(xscl(1)*max(xlims),yscl(1)*max(ylims),zscl(1)*max(zlims(u,:)), ...
            'TUNDRA', ...
            'FontName','Arial', ...
            'FontSize',10);
    end
    
    % Add z-label
    if u == 2
        text(24.5,845,0.25, ...
            'Probability of fire occurring in 30 years (per pixel)', ...
            'FontName','Arial', ...
            'FontSize',8, ...
            'Rotation',90, ...
            'HorizontalAlignment','Center');
    end
    
    % Add x- and y-axis labels
    if u == 1
        text(14,850,-0.02,'T_W_A_R_M (\circC)', ...
            'FontSize',8);
        text(-4,600,-0.012,'P - PET_A_N_N (mm)', ...
            'FontSize',8);
    end
    
    % Format axes
    set(gca, ...
        'FontName','Arial', ...
        'FontSize',8, ...
        'XTick',[2 6 10 14 18], ...
        'YTick',[-500 -250 0 250 500], ...
        'ZTick',[range(zlims(u,:))/nticks(u):range(zlims(u,:))/nticks(u):zlims(u,2)], ...
        'LineWidth',0.72, ...
        'TickDir','in', ...
        'XGrid','on', ...
        'YGrid','on', ...
        'ZGrid','on');
    
    % Rotate view of 3-d surface plot
    view(-137,20);
    
end

% SAVE FIGURE 5
cd(save_dir);
set(gcf,'PaperType','usletter','PaperPositionMode','auto');
print('FIG_5','-dtiff','-r600');
print('FIG_5','-djpeg','-r600');
saveas(gcf,'FIG_5.fig');