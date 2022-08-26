% FIG_4.m
% Matlab version: R2012a
%
% This script creates Figure 4 in Young et al. 2016.
%
% Fig. 4: Partial dependence plots illustrating the relationships between
% the most important explanatory variables and the 30-yr predicted
% probability of fire occurrence. Rows separate different models, with the
% Alaska (AK), boreal forest (BOREAL), and tundra (TUNDRA) models displayed
% from top to bottom. The solid black lines represent the median predicted
% probability of fire occurrence, and the dashed lines represent the
% interquartile range from 100 boosted regression tree models. Probability
% values (y axis) are presented only for the range of climate conditions
% (x axis) observed from 1950-2009. A lowess function (span = 0.1) was used
% to smooth the plotted predicted median and interquartile lines. Vertical
% lines highlight thresholds, identified as the mean breakpoint from the
% segmented regression analysis. As a reference, lighter (darker) colored
% histograms represent the historical distribution of each climate variable
% among unburned (burned) pixels from 1950 to 2009. Histograms heights were
% scaled individually and are not associated with y-axis values.
%
% FILE REQUIREMENTS:
%    Available in '\GBM_Modelling_Analysis_Scripts\AncillaryData\' folder 
%    in 'Figures' nested dataset.
%    ----------------------------------------------------------------------
%    - 'fire_occ_1950_2009.tif' - fire presence/absence information per
%      pixel form 1950 - 2009.
%    - 'TempWarm_1950_2009.tif' - 1950-2009 temperature of the warmest month
%    - 'AnnDEF_1950_2009.tif' - 1950-2009 annual moisture availability
%    - 'masks.mat' - spatial masks for each spatial domain
%
%    Available in '\AK_FINAL_RESULTS\' folder (for example, also BOREAL and
%    TUNDRA).
%    ----------------------------------------------------------------------
%    - 'partDep_TempWarm.csv' - TempWarm partial dependence results from BRTs
%    - 'partDep_AnnDEF.csv' - AnnDEF partial dependence results from BRTs
%
% DEPENDENCIES:
%    * From curve fitting toolbox:
%        - smooth.m
%    * From statistics toolbox:
%        - prctile.m
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

save_dir = [wdir,'Figures\'];

% Load data to create histograms in Fig. 4
cd([wdir,'GBM_Modelling_Analysis_Scripts\AncillaryData\']);
fire_map = geotiffread('fire_occ_1950_2009.tif');

TempWarm_map = ...
    geotiffread('TempWarm_1950_2009.tif');

AnnDEF_map = ...
    geotiffread('AnnDEF_1950_2009.tif');

load('masks.mat');
akmask = masks.akmask;
boreal = masks.boreal;
tundra = masks.tundra;

% Create vectors to plot histograms of fire occurrence and fire absence for
% for each explanatory variable

% Indices of pixels in each spatial domain
idx     = find(akmask(:) == 1);
idx_bor = find(boreal(:) == 1);
idx_tun = find(tundra(:) < 1);

ak_tjja_all   = TempWarm_map(idx);
bor_tjja_all  = TempWarm_map(idx_bor);
tun_tjja_all  = TempWarm_map(idx_tun);

fire_all      = fire_map(idx);
fire_zero     = NaN(size(fire_all));
fire_one      = NaN(size(fire_all));

fire_one(fire_all==1)  = 1;
fire_zero(fire_all==0) = 1;

ak_twrm       = TempWarm_map(akmask==1 & fire_map==0);
ak_twrm_fire  = TempWarm_map(akmask==1 & fire_map==1);
bor_twrm      = TempWarm_map(boreal==1  & fire_map==0);
bor_twrm_fire = TempWarm_map(boreal==1  & fire_map==1);
tun_twrm      = TempWarm_map(tundra==1   & fire_map==0);
tun_twrm_fire = TempWarm_map(tundra==1   & fire_map==1);

ak_adef       = AnnDEF_map(akmask==1 & fire_map==0);
ak_adef_fire  = AnnDEF_map(akmask==1 & fire_map==1);
bor_adef      = AnnDEF_map(boreal==1 & fire_map==0);
bor_adef_fire = AnnDEF_map(boreal==1 & fire_map==1);
tun_adef      = AnnDEF_map(tundra==1 & fire_map==0);
tun_adef_fire = AnnDEF_map(tundra==1 & fire_map==1);

% Parameters to create figure window and panels within figure window
fig_lim_r = 28;
fig_lim_l = 2;
fig_ychg  = 15.75;
fig_xchg  = 16.9;
dxdy = 4.6;
xstart = 3.35;
ystart = 1.5;
btwn_figs_y = 0.2;
btwn_figs_x = 1.0;
nTicks_twrm = 8;
nTicks_adef = 6;

% Initialize figure 4
figure(4);
set(gcf,'Units','Centimeters', ...
    'Position',[fig_lim_r fig_lim_l fig_xchg fig_ychg], ...
    'Color',[1.0 1.0 1.0]);

% -------------------------------------------------------------------------
% ----------------- Partial dependence results for AK ---------------------
% -------------------------------------------------------------------------

% ------------------------- TempWarm --------------------------------------

% Load data from '\AK_FINAL_RESULTS\' folder in 'Generalized Boosted 
% Models and Analysis Scripts' nested dataset. Directory path needs to be
% completed below.
cd([wdir,'AK_FINAL_RESULTS\']);
TempWarm_ak = importdata('partDep_TempWarm.csv');
TempWarm_ak = TempWarm_ak.data;

% Set axes for plotting
axes('units','centimeters','position', ...
    [xstart ystart+(2*dxdy)+(2*btwn_figs_y) dxdy dxdy]);
ylims = [0 0.24];
ylims = [min(ylims):(abs(max(ylims)-min(ylims))/(6)):max(ylims)];
xlims = [4 20];
xlims = [min(xlims):(abs(max(xlims)-min(xlims))/(nTicks_twrm)):max(xlims)];

% Include only those data in the Partial dependence plots under which
% historical climate conditions existed during the 20th-century
include = find(TempWarm_ak(:,1) >= min(ak_twrm) & ...
    TempWarm_ak(:,1) <= max(ak_twrm));

TempWarm_ak = TempWarm_ak(include,:);

Z_ak  = median(TempWarm_ak(:,2:end),2); % Calculate median from 100 BRTs
span = 0.1; % spanning width for lowess smoother

% get interquartile values
[Z_pc] = prctile(TempWarm_ak(:,2:end)',[25 75]);
Z_up_pc = Z_pc(2,:)';
Z_lo_pc = Z_pc(1,:)';

% use lowess smoother on partial dependence data
Zsmooth_ak = smooth(Z_ak,span,'lowess');
Zsmooth_up_pc = smooth(Z_up_pc,span,'lowess');
Zsmooth_lo_pc = smooth(Z_lo_pc,span,'lowess');

hold on;

% Plot histogram results
nBins = [2:0.5:20];
[n,xout] = hist(ak_twrm,nBins);
bh = bar(xout',(n/(size(ak_twrm_fire,1)+size(ak_twrm,1)))'*2,'hist');
set(bh(1),'FaceColor',[.7 .7 .7],'EdgeColor','none');

[n,xout] = hist(ak_twrm_fire,nBins);
bh = bar(xout',(n/(size(ak_twrm_fire,1)+size(ak_twrm,1)))'*2,'hist');
set(bh(1),'FaceColor',[.3 .3 .3],'EdgeColor','none');

% Plot median and interquartile predicted probability of fire occurrence
% from the distribution of the 100 BRTs
plot(TempWarm_ak(:,1),Zsmooth_up_pc,'-.','Color','k','LineWidth',1.44);
plot(TempWarm_ak(:,1),Zsmooth_lo_pc,'-.','Color','k','LineWidth',1.44);
plot(TempWarm_ak(:,1),Zsmooth_ak,'-', ...
    'LineWidth',1.44,'Color','k'); hold on;

% Set axis limits
ylim([min(ylims) max(ylims)]);
xlim([min(xlims) max(xlims)]);

axis square;
box on;

% Add text to figure panel
text(4.6,0.9*max(ylims),'(a)', ...
    'FontName','Arial', ...
    'FontWeight','Bold');

text(0.38*range(xlims)+min(xlims),0.9*max(ylims),'AK', ....
    'FontName','Arial', ...
    'FontSize',9,...
    'HorizontalAlignment','Center');

% Format axes
set(gca,'XTick',xlims,...
    'XTickLabel',{''},...
    'YTick',ylims(2:(end-1)), ...
    'FontName','Arial',...
    'FontSize',8, ...
    'LineWidth',0.72);
XTick_data = [(xlims(2:(end-1)));
    (xlims(2:(end-1)))];
YTick_data = [zeros(1,length(xlims(2:(end-1))));
    ((0.01*max(ylims))')*ones(1,length(xlims(2:(end-1))))];

% Add tick lines to x-axis - unable to see in current plot due to
% histograms.
line(XTick_data,YTick_data,'Color','k','LineWidth',0.72)

% Add line indicating mean estimated climatic threshold(s)
line([13.36 13.36],[0,.5],'Color','k','LineStyle','-',...
    'LineWidth',1.44);

% ------------------------- AnnDEF ----------------------------------------

% Load data
AnnDef_ak = importdata('partDep_AnnDEF.csv');
AnnDef_ak = AnnDef_ak.data;

% Set axes for plotting
axes('units','centimeters','position', ...
    [xstart+dxdy+btwn_figs_x ystart+(2*dxdy)+(2*btwn_figs_y) dxdy dxdy]);

ylims = [0 0.24];
ylims = [min(ylims):(abs(max(ylims)-min(ylims))/(6)):max(ylims)];
xlims = [-450 450];
xlims = [min(xlims):(abs(max(xlims)-min(xlims))/(nTicks_adef)):max(xlims)];

% Include only those data in the Partial dependence plots under which
% historical climate conditions existed during the 20th-century
include = find(AnnDef_ak(:,1) >= min(ak_adef) & ...
    AnnDef_ak(:,1) <= max(ak_adef));

AnnDef_ak = AnnDef_ak(include,:);

Z_ak  = median(AnnDef_ak(:,2:end),2); % Calculate median from 100 BRTs
span = 0.1; % spanning width for lowess smoother

% get interquartile values
[Z_pc] = prctile(AnnDef_ak(:,2:end)',[25 75]);
Z_up_pc = Z_pc(2,:)';
Z_lo_pc = Z_pc(1,:)';

% use lowess smoother on partial dependence data
Zsmooth_ak = smooth(Z_ak,span,'lowess');
Zsmooth_up_pc = smooth(Z_up_pc,span,'lowess');
Zsmooth_lo_pc = smooth(Z_lo_pc,span,'lowess');

hold on;

% Plot histogram results
nBins = [-500:25:10000];
[n,xout] = hist(ak_adef,nBins);
bh = bar(xout',(n/(size(ak_adef_fire,1)+size(ak_adef,1)))'*3.5,'hist');
set(bh(1),'FaceColor',[.7 .7 .7],'EdgeColor','none');

[n,xout] = hist(ak_adef_fire,nBins);
bh = bar(xout',(n/(size(ak_adef_fire,1)+size(ak_adef,1)))'*3.5,'hist');
set(bh(1),'FaceColor',[.3 .3 .3],'EdgeColor','none');

% Plot median and interquartile predicted probability of fire occurrence
% from the distribution of the 100 BRTs
plot(AnnDef_ak(:,1),Zsmooth_up_pc,'-.','Color','k','LineWidth',1.44);
plot(AnnDef_ak(:,1),Zsmooth_lo_pc,'-.','Color','k','LineWidth',1.44);
plot(AnnDef_ak(:,1),Zsmooth_ak,'-', ...
    'LineWidth',1.44,'Color','k');

% Set axis limits
xlim([min(xlims) max(xlims)]);
ylim([min(ylims) max(ylims)]);

% Add text to figure panel
text(0.7*max(xlims),0.9*max(ylims),'(b)','FontName','Arial', ...
    'FontSize',11,'FontWeight','Bold');
text(0.5*range(xlims)+min(xlims),0.9*max(ylims),'AK', ....
    'FontName','Arial', ...
    'FontSize',9,...
    'HorizontalAlignment','Center');

% Format axes
axis square;
box on;
set(gca,'XTick',xlims,...
    'XTickLabel',{''},...
    'YTick',ylims(2:(end-1)),...
    'FontName','Arial', ...
    'LineWidth',0.72, ...
    'FontSize',8);
XTick_data = [xlims(2:(end-1));
    xlims(2:(end-1))];
YTick_data = [zeros(1,length(xlims(2:(end-1))));
    ((0.01*max(ylims))')*ones(1,length(xlims(2:(end-1))))];

% Add tick lines to x-axis - unable to see in current plot due to
% histograms.
line(XTick_data,YTick_data,'Color','k','LineWidth',0.72)

% Add line indicating mean estimated climatic threshold(s)
line([215 215],[0,.24],'Color','k','LineStyle','-',...
    'LineWidth',1.44);

% -------------------------------------------------------------------------
% ------------- Partial dependence results for BOREAL models --------------
% -------------------------------------------------------------------------

% ------------------------- TempWarm --------------------------------------

% Load data from '\BOREAL_FINAL_RESULTS\' folder in 'Generalized Boosted 
% Models and Analysis Scripts' nested dataset. Directory path needs to be
% completed below.
cd([wdir,'BOREAL_FINAL_RESULTS\']);
TempWarm_bor = importdata('partDep_TempWarm.csv');
TempWarm_bor = TempWarm_bor.data;

% Set axes for plotting
axes('units','centimeters','position', ...
    [xstart ystart+(1*dxdy)+(1*btwn_figs_y) dxdy dxdy]);
ylims = [0 0.30];
ylims = [min(ylims):(abs(max(ylims)-min(ylims))/(6)):max(ylims)];
xlims = [4 20];
xlims = [min(xlims):(abs(max(xlims)-min(xlims))/(nTicks_twrm)):max(xlims)];

% Include only those data in the Partial dependence plots under which
% historical climate conditions existed during the 20th-century
include = find(TempWarm_bor(:,1) >= min(bor_twrm) & ...
    TempWarm_bor(:,1) <= max(bor_twrm));

TempWarm_bor = TempWarm_bor(include,:);

Z_bor  = median(TempWarm_bor(:,2:end),2); % Calculate median from 100 BRTs

% get interquartile values
[Z_pc] = prctile(TempWarm_bor(:,2:end)',[25 75]);
Z_up_pc = Z_pc(2,:)';
Z_lo_pc = Z_pc(1,:)';

% use lowess smoother on partial dependence data
Zsmooth_bor = smooth(Z_bor,span,'lowess');
Zsmooth_up_pc = smooth(Z_up_pc,span,'lowess');
Zsmooth_lo_pc = smooth(Z_lo_pc,span,'lowess');

% Plot histogram results
nBins = [2:0.5:20];
[n,xout] = hist(bor_twrm,nBins);
bh = bar(xout',(n/(size(bor_twrm_fire,1)+size(bor_twrm,1)))'*3,'hist');
set(bh(1),'FaceColor',[.7 .7 .7],'EdgeColor','none');

hold on;

% Plot histogram results
[n,xout] = hist(bor_twrm_fire,nBins);
bh = bar(xout',(n/(size(bor_twrm_fire,1)+size(bor_twrm,1)))'*3,'hist');
set(bh(1),'FaceColor',[.3 .3 .3],'EdgeColor','none');

% Plot median and interquartile predicted probability of fire occurrence
% from the distribution of the 100 BRTs
plot(TempWarm_bor(:,1),Zsmooth_up_pc,'-.','Color','k','LineWidth',1.44);
plot(TempWarm_bor(:,1),Zsmooth_lo_pc,'-.','Color','k','LineWidth',1.44);
plot(TempWarm_bor(:,1),Zsmooth_bor,'-', ...
    'LineWidth',1.44, ...
    'Color','k');

% Set axis limits
xlim([min(xlims) max(xlims)]);
ylim([min(ylims) max(ylims)]);

box on;
axis square;

% Add text to figure panel
text(4.6,0.9*max(ylims),'(c)', ...
    'FontName','Arial', ...
    'FontSize',11, ...
    'FontWeight','Bold');

text(0.38*range(xlims)+min(xlims),0.9*max(ylims),'BOREAL', ....
    'FontName','Arial', ...
    'FontSize',9,...
    'HorizontalAlignment','Center');

% Format axes
set(gca,'XTick',xlims,...
    'XTickLabel',{''}, ...
    'YTick',ylims(2:(end-1)),...
    'FontName','Arial', ...
    'FontSize',8, ...
    'LineWidth',0.72);
XTick_data = [xlims(2:(end-1));
    xlims(2:(end-1))];
YTick_data = [zeros(1,length(xlims(2:(end-1))));
    ((0.01*max(ylims))')*ones(1,length(xlims(2:(end-1))))];

% Add tick lines to x-axis - unable to see in current plot due to
% histograms.
line(XTick_data,YTick_data,'Color','k','LineWidth',0.72)

% Add line indicating mean estimated climatic threshold(s)
line([13.47 13.47],[0,.5],'Color','k','LineStyle','-',...
    'LineWidth',1.44);

% Add y-axis label
ylabel('Probability of fire cccurring in 30 years (per pixel)', ...
    'FontName','Arial', ...
    'FontSize',10);
% --------------------------- AnnDEF --------------------------------------

% Load data
AnnDef_bor = importdata('partDep_AnnDEF.csv');
AnnDef_bor = AnnDef_bor.data;

% Set axes for plotting
axes('units','centimeters','position', ...
    [xstart+dxdy+btwn_figs_x ystart+(1*dxdy)+(1*btwn_figs_y) dxdy dxdy]);
ylims = [0 0.30];
ylims = [min(ylims):(abs(max(ylims)-min(ylims))/(6)):max(ylims)];
xlims = [-450 450];
xlims = [min(xlims):(abs(max(xlims)-min(xlims))/(nTicks_adef)):max(xlims)];

% Include only those data in the Partial dependence plots under which
% historical climate conditions existed during the 20th-century
include = find(AnnDef_bor(:,1) >= min(bor_adef) & ...
    AnnDef_bor(:,1) <= max(bor_adef));

AnnDef_bor = AnnDef_bor(include,:);

Z_bor  = median(AnnDef_bor(:,2:end),2); % Calculate median from 100 BRTs

% get interquartile values
[Z_pc] = prctile(AnnDef_bor(:,2:end)',[25 75]);
Z_up_pc = Z_pc(2,:)';
Z_lo_pc = Z_pc(1,:)';

% use lowess smoother on partial dependence data
Zsmooth_bor   = smooth(Z_bor,span,'lowess');
Zsmooth_up_pc = smooth(Z_up_pc,span,'lowess');
Zsmooth_lo_pc = smooth(Z_lo_pc,span,'lowess');

% Plot histogram results
nBins = [-500:25:10000];
[n,xout] = hist(bor_adef,nBins);
bh = bar(xout',(n/(size(bor_adef_fire,1)+size(bor_adef,1)))'*4.5,'hist');
set(bh(1),'FaceColor',[.7 .7 .7],'EdgeColor','none');

hold on;

[n,xout] = hist(bor_adef_fire,nBins);
bh = bar(xout',(n/(size(bor_adef_fire,1)+size(bor_adef,1)))'*4.5,'hist');
set(bh(1),'FaceColor',[.3 .3 .3],'EdgeColor','none');

% Plot median and interquartile predicted probability of fire occurrence
% from the distribution of the 100 BRTs
plot(AnnDef_bor(:,1),Zsmooth_up_pc,'-.','Color','k','LineWidth',1.44);
plot(AnnDef_bor(:,1),Zsmooth_lo_pc,'-.','Color','k','LineWidth',1.44);
plot(AnnDef_bor(:,1),Zsmooth_bor, ...
    '-','LineWidth',1.44, ...
    'Color','k');

% Set axis limits
ylim([min(ylims) max(ylims)]);
xlim([min(xlims) max(xlims)]);

axis square;
box on;

% Add text to figure panel
text(0.70*max(xlims),0.90*max(ylims),'(d)', ...
    'FontName','Arial', ...
    'FontSize',11, ...
    'FontWeight','Bold');

text(0.5*range(xlims)+min(xlims),0.9*max(ylims),'BOREAL', ....
    'FontName','Arial', ...
    'FontSize',9,...
    'HorizontalAlignment','Center');

% Format axes
set(gca,'XTick',xlims,...
    'XTickLabel',{''}, ...
    'YTick',ylims(2:(end-1)),...
    'FontName','Arial', ...
    'FontSize',8, ...
    'LineWidth',0.72);
XTick_data = [xlims(2:(end-1));
    xlims(2:(end-1))];
YTick_data = [zeros(1,length(xlims(2:(end-1))));
    ((0.01*max(ylims))')*ones(1,length(xlims(2:(end-1))))];

% Add tick lines to x-axis - unable to see in current plot due to
% histograms.
line(XTick_data,YTick_data,'Color','k','LineWidth',0.72)

% Add line indicating mean estimated climatic threshold(s)
line([151 151],[0,0.5],'Color','k','LineStyle','-',...
    'LineWidth',1.44);

% -------------------------------------------------------------------------
% ------------- Partial dependence results for TUNDRA models --------------
% -------------------------------------------------------------------------

% ------------------------- TempWarm --------------------------------------

% Load data from '\TUNDRA_FINAL_RESULTS\' folder in 'Generalized Boosted 
% Models and Analysis Scripts' nested dataset. Directory path needs to be
% completed below.
cd([wdir,'TUNDRA_FINAL_RESULTS\']);
TempWarm_tun = importdata('partDep_TempWarm.csv');
TempWarm_tun = TempWarm_tun.data;

% Set axes for plotting
axes('units','centimeters','position',[xstart ystart dxdy dxdy]);
ylims = [0 0.12];
ylims = [min(ylims):(abs(max(ylims)-min(ylims))/(6)):max(ylims)];
xlims = [4 20];
xlims = [min(xlims):(abs(max(xlims)-min(xlims))/(nTicks_twrm)):max(xlims)];

% Include only those data in the Partial dependence plots under which
% historical climate conditions existed during the 20th-century
include = find(TempWarm_tun(:,1) >= min(tun_twrm) & ...
    TempWarm_tun(:,1) <= max(tun_twrm));

TempWarm_tun = TempWarm_tun(include,:);

Z_tun  = median(TempWarm_tun(:,2:end),2); % Calculate median from 100 BRTs

% get interquartile values
[Z_pc] = prctile(TempWarm_tun(:,2:end)',[25 75]);
Z_up_pc = Z_pc(2,:)';
Z_lo_pc = Z_pc(1,:)';

% use lowess smoother on partial dependence data
Zsmooth_tun = smooth(Z_tun,span,'lowess');
Zsmooth_up_pc = smooth(Z_up_pc,span,'lowess');
Zsmooth_lo_pc = smooth(Z_lo_pc,span,'lowess');

% Plot histogram results
nBins = [2:0.5:20];
[n,xout] = hist(tun_twrm,nBins);
bh = bar(xout',(n/(size(tun_twrm_fire,1)+size(tun_twrm,1)))'*0.6,'hist');
set(bh(1),'FaceColor',[.7 .7 .7],'EdgeColor','none');

hold on;

[n,xout] = hist(tun_twrm_fire,nBins);
bh = bar(xout',(n/(size(tun_twrm_fire,1)+size(tun_twrm,1)))'*0.6,'hist');
set(bh(1),'FaceColor',[.3 .3 .3],'EdgeColor','none');

% Plot median and interquartile predicted probability of fire occurrence
% from the distribution of the 100 BRTs
plot(TempWarm_tun(:,1),Zsmooth_up_pc,'-.','Color','k', ...
    'LineWidth',1.44);
plot(TempWarm_tun(:,1),Zsmooth_lo_pc,'-.','Color','k', ...
    'LineWidth',1.44);

plot(TempWarm_tun(:,1),Zsmooth_tun,'-', ...
    'LineWidth',2, ...
    'Color','k', ...
    'LineWidth',1.44);

% Set axis limits
ylim([min(ylims) max(ylims)])
xlim([min(xlims) max(xlims)])

axis square;
box on;

% Add text to figure panel
text(4.6,0.9*max(ylims),'(e)', ...
    'FontName','Arial', ...
    'FontSize',11, ...
    'FontWeight','Bold');

text(0.38*range(xlims)+min(xlims),0.9*max(ylims),'TUNDRA', ....
    'FontName','Arial', ...
    'FontSize',9,...
    'HorizontalAlignment','Center');

% Format axes
set(gca,'XTick',xlims,...
    'YTick',ylims(2:(end-1)),...
    'FontName','Arial', ...
    'FontSize',8, ...
    'LineWidth',0.72);
XTick_data = [xlims(2:(end-1));
    xlims(2:(end-1))];
YTick_data = [zeros(1,length(xlims(2:(end-1))));
    ((0.01*max(ylims))')*ones(1,length(xlims(2:(end-1))))];

% Add tick lines to x-axis - unable to see in current plot due to
% histograms.
line(XTick_data,YTick_data,'Color','k','LineWidth',0.72)

% Add line indicating mean estimated climatic threshold(s)
line([13.6 13.6],[0,0.5],'Color','k','LineStyle','-',...
    'LineWidth',1.44);

% Add x-axis label
xlabel({'Mean temp. of the', 'warmest month ( \circC)'}, ...
    'FontSize',10);

% ------------------------- AnnDEF ----------------------------------------

% Load data
AnnDef_tun = importdata('partDep_AnnDEF.csv');
AnnDef_tun = AnnDef_tun.data;

% Set axes for plotting
axes('units','centimeters','position', ...
    [xstart+dxdy+btwn_figs_x ystart dxdy dxdy]);
ylims = [0 0.018];
ylims = [min(ylims):(abs(max(ylims)-min(ylims))/(6)):max(ylims)];
xlims = [-450 450];
xlims = [min(xlims):(abs(max(xlims)-min(xlims))/(nTicks_adef)):max(xlims)];

% Include only those data in the Partial dependence plots under which
% historical climate conditions existed during the 20th-century
include = find(AnnDef_tun(:,1) >= min(tun_adef) & ...
    AnnDef_tun(:,1) <= max(tun_adef));

AnnDef_tun = AnnDef_tun(include,:);

Z_tun  = median(AnnDef_tun(:,2:end),2); % Calculate median from 100 BRTs

% get interquartile values
[Z_pc] = prctile(AnnDef_tun(:,2:end)',[25 75]);
Z_up_pc = Z_pc(2,:)';
Z_lo_pc = Z_pc(1,:)';

% use lowess smoother on partial dependence data
Zsmooth_tun   = smooth(Z_tun,span,'lowess');
Zsmooth_up_pc = smooth(Z_up_pc,span,'lowess');
Zsmooth_lo_pc = smooth(Z_lo_pc,span,'lowess');

% Plot histogram results
nBins = [-500:25:10000];
[n,xout] = hist(tun_adef,nBins);
bh = bar(xout',(n/(size(tun_adef_fire,1)+size(tun_adef,1)))'*.2,'hist');
set(bh(1),'FaceColor',[.7 .7 .7],'EdgeColor','none');

hold on;

[n,xout] = hist(tun_adef_fire,nBins);
bh = bar(xout',(n/(size(tun_adef_fire,1)+size(tun_adef,1)))'*.2,'hist');
set(bh(1),'FaceColor',[.3 .3 .3],'EdgeColor','none');

% Plot median and interquartile predicted probability of fire occurrence
% from the distribution of the 100 BRTs
plot(AnnDef_tun(:,1),Zsmooth_up_pc,'-.','Color','k','LineWidth',1.44);
plot(AnnDef_tun(:,1),Zsmooth_lo_pc,'-.','Color','k','LineWidth',1.44);
plot(AnnDef_tun((1:(end-2)),1),Zsmooth_tun(1:(end-2)),'-', ...
    'LineWidth',1.44,...
    'Color','k');

% Set axis limits
xlim([min(xlims) max(xlims)]);
ylim([min(ylims) max(ylims)]);


% Add text to figure panel
text(0.70*max(xlims),0.9*max(ylims),'(f)', ...
    'FontName','Arial', ...
    'FontWeight','Bold', ...
    'FontSize',11);

text(0.5*range(xlims)+min(xlims),0.9*max(ylims),'TUNDRA', ....
    'FontName','Arial', ...
    'FontSize',9,...
    'HorizontalAlignment','Center');

% Format axes
box on;
axis square;
set(gca,'XTick',xlims(2:(end-1)), ...
    'YTick',ylims(2:(end-1)), ...
    'FontName','Arial', ...
    'FontSize',8, ...
    'LineWidth',0.72);
XTick_data = [xlims(2:(end-1));
    xlims(2:(end-1))];
YTick_data = [zeros(1,length(xlims(2:(end-1))));
    ((0.01*max(ylims))')*ones(1,length(xlims(2:(end-1))))];

% Add tick lines to x-axis - unable to see in current plot due to
% histograms.
line(XTick_data,YTick_data,'Color','k','LineWidth',0.72)

% Add line(s) indicating mean estimated climatic threshold(s)
line([154 154],[0,.5],'Color','k','LineStyle','-',...
    'LineWidth',1.44);
line([-207 -207],[0,.5],'Color','k','LineStyle','-',...
    'LineWidth',1.44);

% Add x-axis label
xlabel({'Mean total annual','P - PET (mm)'},...
    'FontSize',10,'FontName','Arial');

% SAVE FIGURE AS .FIG AND .TIF FILES AT 600 DPI RESOLUTION
cd(save_dir);
set(gcf,'PaperType','usletter','PaperPositionMode','auto');
print('FIG_4','-dtiff','-r600');
print('FIG_4','-djpeg','-r600');
saveas(gcf,'FIG_4.fig');