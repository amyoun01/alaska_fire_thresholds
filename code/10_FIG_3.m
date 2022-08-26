% FIG_3.m
% Matlab version: R2012a
%
% This script creates Figure 3 in Young et al. 2016. 
% 
% Fig. 3: Relative influence of explanatory variables for the Alaska (AK), 
% boreal forest (BOREAL), and tundra (TUNDRA) models. Bar heights represent 
% the sample means and error bars represent ± 1 standard deviation from 100 
% boosted regression tree models. For the BOREAL model, the relative 
% influence of vegetation (Veg) is 0 by default, as the BOREAL vegetation 
% model has only one class (indicated by the black diamond). 
%
% FILE REQUIREMENTS:
%
%    Available in '\AK_FINAL_RESULTS\', '\BOREAL_FINAL_RESULTS\', and 
%    '\TUNDRA_FINAL_RESULTS\' folders in 'Generalized Boosted 
%    Models and Analysis Scripts' nested dataset.
%    ----------------------------------------------------------------------
%    - 'relInf.csv' - Relative influence results from BRT analysis
%
% DEPENDENCIES:
%    * barwitherr.m - function written by Martina Callaghan and available
%        from the Mathworkds File Exchange:
%        'http://uk.mathworks.com/matlabcentral/fileexchange/30639-bar-chart-with-error-bars/content/barwitherr.m'
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
                % nested datasets are stored. Labeled as [WDIR] in 
                % remainder of script. This may need to be further modified
                % below depending on how directories are organized by user.

% Directories for results from each of the three models
modeltype = {'1_AK_FINAL_RESULTS', ...
             '2_BOREAL_FINAL_RESULTS', ...
             '3_TUNDRA_FINAL_RESULTS'};
% Directory to save figure to
save_dir = [wdir,'Figures\'];

n = 4; % NUMBER OF EXPLANATORY VARIABLES USED IN BRT MODELS

mRel = NaN*ones(n,length(modeltype)); % ALLOCATED SPACE TO STORE MEAN
                                      % RELATIVE INFLUENCE VALUES FROM THE
                                      % 100 BRTS FOR EACH MODEL (AK,
                                      % BOREAL,TUNDRA)
                                      
sRel = NaN*ones(n,length(modeltype)); % ALLOCATED SPACE TO STORE SAMPLE
                                      % STANDARD DEVIATION VALUES OF
                                      % RELATIVE INFLUENCE ESTIMATES FROM 
                                      % THE 100 BRTS FOR EACH MODEL (AK,
                                      % BOREAL,TUNDRA)

for m = 1:length(modeltype) % FOR THE THREE MODELS
    
    % CREATE DIRECTORY NAME TO WHERE FILES ARE STORED.
    chgdir    = [wdir,char(modeltype(m))];
    cd(chgdir);
    
    % IMPORT RELATIVE INLFUENCE DATA
    relInf  = importdata('relInf.csv');
    
    mRel(:,m) = mean(relInf.data); % mean relative influence for each 
                                   % explanatory variable
    sRel(:,m) = std(relInf.data); % SD relative influence for each 
                                   % explanatory variable
    
end

% Initialize figure
figure(3); clf; 
set(gcf, ...
    'Units','Centimeters', ...
    'Position',[27 9 8 8.8]);

% Labels for bar graph
labelsB = {'P - PET_A_N_N','T_W_A_R_M','TR','Veg'};

% order to plot bar graphs in
idx = [2 1 3 4];

% create bar graph
[h hE] = barwitherr(sRel(idx,:),mRel(idx,:),'LineWidth',1);

set(hE(1),'LineWidth',1.44);
set(hE(2),'LineWidth',1.44);
set(hE(3),'LineWidth',1.44);

set(gca,'XTickLabel',{}, ...
        'YTickLabel',{'10','20','30','40','50','60','70'}, ...
        'FontName','Arial', ...
        'FontSize',10, ...
        'LineWidth',0.72, ...
        'Units','Centimeters', ...
        'Position',[1.2 2 6.7 6.7], ...
        'XTick',[], ...
        'YTick',[10:10:70]);
axis square;

% SET X-AXIS AND Y-AXIS LIMITS
xlim([0 4.5]);
ylim([0 80]); 

% ADD YLABEL TO PLOT
ylabel('Relative influence (%)', ...
       'FontSize',10, ...
       'FontName','Arial', ...
       'FontWeight','Normal');

% COLORMAP FOR BAR PLOTS DISPLAYING RELATIVE INFLUENCE DATA
cMap = [0.50 0.50 0.50;  % AK
        0.75 0.75 0.75;  % BOREAL
        1.00 1.00 1.00]; % TUNDRA
    
colormap(cMap); % APPLY COLORMAP

% ADD LEGEND
lg = legend('AK','BOREAL','TUNDRA');
legend('boxoff'); % REMOVE BOX FROM LEGEND

pos = get(gca,'Position');

for i = 1:n
   text(i,pos(2)-5,labelsB(idx(i)), ...
       'FontWeight','Normal', ...
       'FontSize',10, ...
       'FontName','Arial', ...
       'Rotation',300) 
end

hold on;

% Add diamond to veg type relative influence for borel forest
plot(4,6,'kd', ...
     'MarkerFaceColor','k');

% SAVE FIGURE AS .FIG AND .TIF FILES AT 600 DPI RESOLUTION
cd(save_dir);
set(gcf,'PaperType','usletter','PaperPositionMode','auto');
print('FIG_3','-dtiff','-r600');
print('FIG_3','-djpeg','-r600');
saveas(gcf,'FIG_3.fig');