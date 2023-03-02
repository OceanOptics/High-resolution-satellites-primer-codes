%------------------------------------------------------------------------ %
%                                                                         %
% MATLAB CODE TO PLOT data from SATELLITE SCENES (L5, L7, L8 and L9) from %
% Patos Lagoon ESTUARY                                                    %                                                 
%                                                                         %
% THIS CODE REQUIRES:                                                     %
%                                                                         %
% Inputs:                                                                 %
% 1)'Landsat_matchups.mat' file                                           %
% 2)'Turb_Buoy4.mat'
% 3)'Turb_L5L7L8L9_histogram.mat')                                        %
% 3) L5 scenes (L5_TM_1985_08_17_12_48_44_221082_L2W.nc)                  %
% 4) L7 scenes (L7_ETM_2002_09_09_13_07_03_221082_L2W.nc)                 %
% 5) L8 scenes (L8_OLI_2014_01_21_13_20_15_221082_L1R.nc)                 %
% 6) L9 scenes (L5_TM_1985_08_17_12_48_44_221082_L2W.nc)                  %
% 7) m_map package to plot maps
% 8) lsqbisec.m function from MBARI for type 2 regression
%                                                                         %
%                                                                         %
% developed by:                                                           %
% Juliana Tavora (j.tavora@utwente.nl)                                    %
% University of Twente                                                    %
% version November 2022                                                   %
%                                                                         %
%------------------------------------------------------------------------ %

clear all
close all
clc


load('Landsat_matchups.mat')

%-------------------------------------------------------------------------%
%  step 1 -  matchups with in-situ data (Buoy RS1 and RS4) from SIMCOSTA (BRASIL)

% determine in-situ stations in matchups to regioanlly calibrate and plot %

% here in hands with both in-situ data and satellite sampled data, we used
% MATLABs curve-fitting app tool to find the curve that best approximates
% the response relationship between the turbidity data. then we manually
% collect the coeficients and metrics for further analises. 
% 
% for this exercise, a power-law curve ( y = ax^b ) was used:
%       a = 0.0993
%       b = 1.354
%       x = turbidity derived from satellite
%       y = turbidity from two buoys (in-situ data)
%
%-------------------------------------------------------------------------%

stations = unique(big_table_matchup_l8l9.boias);


[xData, yData] = prepareCurveData(big_table_matchup_l8l9.Avg_Turb,big_table_matchup_l8l9.Turb_Nechad16_865_median);
idx = ismember(big_table_matchup_l8l9.boias,stations(2));

T_sat = yData(~idx); T_insitu = xData(~idx);
T_corrected = 0.0993.*T_sat.^(1.354); 

%--------------------------- plot match-up  ------------------------------%

stations = unique(big_table_matchup_l8l9.boias);
symb={'o','>', '<', 's'};


figure(1)
subaxis(5,5,[1:2 5:6],'Spacing',0.01,'MR',0.1);

for i=[1,3]%:length(stations)
    
    idx = ismember(big_table_matchup_l8l9.boias,stations(i));
    hold all
    
    errorbar(((big_table_matchup_l8l9.Avg_Turb(idx))), ...
        (0.0993.*(big_table_matchup_l8l9.Turb_Nechad16_865_median(idx)).^(1.354)), ...
        (big_table_matchup_l8l9.Turb_Nechad16_865_std(idx)),'k', 'LineStyle','none'); %percentual
    
    dataplot(i) = scatter((big_table_matchup_l8l9.Avg_Turb(idx)), ...
        (0.0993.*(big_table_matchup_l8l9.Turb_Nechad16_865_median(idx)).^(1.354)),50,...
        symb{i},'MarkerFaceColor', 'w' ,'MarkerEdgeColor', 'k' ,'LineWidth',1);
    
    %metrics
    N(i-0)                  = size(~isnan(big_table_matchup_l8l9.Avg_Turb(idx)),1);
    [r_kendall(i-0),p(i-0)] = corr(big_table_matchup_l8l9.Avg_Turb(idx),(0.0993.*(big_table_matchup_l8l9.Turb_Nechad16_865_median(idx)).^(1.354)),'Type','Kendall','Rows','pairwise');
    MAPD(i-0)                = nansum( abs ((0.0993.*(big_table_matchup_l8l9.Turb_Nechad16_865_median(idx)).^(1.354)) - big_table_matchup_l8l9.Avg_Turb(idx)  )./...
                                        big_table_matchup_l8l9.Avg_Turb(idx)  ) ./ N(i-0);
    RMSD(i-0)               = sqrt(nansum(((0.0993.*(big_table_matchup_l8l9.Turb_Nechad16_865_median(idx)).^(1.354))  - big_table_matchup_l8l9.Avg_Turb(idx)).^2 )/ N(i-0));
    
end



%  metrics
idx = ismember(big_table_matchup_l8l9.boias,stations(2));
N(4)                  = size(~isnan(big_table_matchup_l8l9.Avg_Turb(~idx)),1);
[r_kendall(4),p(4)]   = corr(big_table_matchup_l8l9.Avg_Turb(~idx),T_corrected,'Type','Kendall','Rows','pairwise');
MAPD(4)               = nansum((abs(T_corrected - big_table_matchup_l8l9.Avg_Turb(~idx)))./...
    ((big_table_matchup_l8l9.Avg_Turb(~idx)))) ./ N(4);
RMSD(4)               = sqrt(nansum((T_corrected  - big_table_matchup_l8l9.Avg_Turb(~idx)).^2 )/ N(4));


%  type 2 regression
[m,b,~,~,~]  = lsqbisec(big_table_matchup_l8l9.Avg_Turb(~idx),(0.0993.*(big_table_matchup_l8l9.Turb_Nechad16_865_median(~idx)).^(1.354)));
latlot       = linspace(0.001,400,500);
lonlot       = m.*latlot + b;

hold on
oneone = plot(([0.001:400]),([0.001:400]),':k');
reg= plot(latlot,lonlot,'k','LineWidth',0.5);
set(gca,'FontSize',12)
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'linear')
xlim([0.1 120])
ylim([0.1 120])
xlabel('Turbidity in-situ [NTU]','FontSize',12)
ylabel('Turbidity derived [NTU]','FontSize',12)
set(gca,'color',[1 1 1]); %sets grey background

legend([oneone reg dataplot(1)],{'1:1','y = 0.79x + 0.31 (r = 0.45, p = 0.0002, n = 31)'},'EdgeColor','none','Color','none','Location','eastoutside','FontSize',12)% box on

%-------------------------------------------------------------------------%
% step 2 - determine histogram of distribution for regionally calibrated satellite-derived Turb

load('Turb_Buoy4.mat')
load('Turb_L5L7L8L9_histogram.mat')

% plot histogram

figure(1)
subaxis(5,5,[3:4 7:8],'Spacing',0.01,'MR',0.1);

field = histogram(log10(Turb_buoy4.Avg_Turb),50,'Normalization','probability','FaceColor','k');
hold on
set(gca,'FontSize',12)

edge = field.BinEdges;
sat = histogram(log10(0.0993.*(big_table_Landsat_RS4.Turb_NIR_box_mean).^(1.354)),edge,'Normalization','probability','FaceColor',[0.47,0.67,0.19]);
hold on
set(gca,'FontSize',12)
ylim([0 0.1])


leg = legend([field, sat],{'in-situ: 28.1 $\pm$ 44.45,','satellite-derived: 6.26 $\pm$ 44.39'},'EdgeColor','none','interpreter','latex');
title(leg,' Turbidity [FTU]')
% set(gca,'FontSize',12)

xlabel('Turbidity [NTU]','FontSize',12)
ylabel('Relative frequency','FontSize',12)


%-------------------------------------------------------------------------%
% step 3 - plot examples of (regionally calibrated) turbidity maps from each landsat sensor 

pathDir = '/path_to/HighResol_primer/';
list_scenes = dir([pathDir 'L5*L2W.nc']);

for i=2
    
    fileName = char({list_scenes(i).name});
    imageFile = [pathDir fileName];
    SPM = ncread(imageFile,'TUR_Nechad2009_839');  SPM = SPM'; SPM(SPM <=0) = NaN;
    
    %latitude/long from the image
    lat = ncread(imageFile, 'lat'); lat = lat';
    lon = ncread(imageFile, 'lon'); lon = lon';
end

figure(1)
subaxis(5,4,[13 17],'Spacing',0.01,'MR',0.1);
m_proj('mercator','lon',  [-52.3 -51.8], 'lat',[-32.4 -31.7]);
m_pcolor(lon,lat,0.0993.*SPM.^(1.354));
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linewidth', 0.6, 'tickdir', 'none','gridcolor','none', 'backgroundcolor',[0.901960784313726 0.901960784313726 0.901960784313726],'fontsize',12);
m_text(-52.25,-32.35,'L5-TM (2001-05-25)','color','w','fontsize',12);

caxis([0 100])

hold on

list_scenes = dir([pathDir 'L7*L2W.nc']);

for i=4
    fileName = char({list_scenes(i).name});
    imageFile = [pathDir fileName];
    SPM = ncread(imageFile,'TUR_Nechad2009_835');  SPM = SPM'; SPM(SPM <=0) = NaN;
    
    %latitude/long from the image
    lat = ncread(imageFile, 'lat'); lat = lat';
    lon = ncread(imageFile, 'lon'); lon = lon';
end

subaxis(5,4,[14 18],'Spacing',0.01,'MR',0.1);
m_proj('mercator','lon',  [-52.3 -51.8], 'lat',[-32.4 -31.7]);
m_pcolor(lon,lat,0.0993.*SPM.^(1.354));
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linewidth', 0.6, 'tickdir', 'none','gridcolor','none', 'backgroundcolor',[0.901960784313726 0.901960784313726 0.901960784313726],'fontsize',12);
m_text(-52.25,-32.35,'L7-ETM (2002-09-09)','color','w','fontsize',12);

caxis([0 100])

hold on

list_scenes = dir([pathDir 'L8*L2W.nc']);

for i=1
    
    fileName = char({list_scenes(i).name});
    imageFile = [pathDir fileName];
    SPM = ncread(imageFile,'TUR_Nechad2009_865');  SPM = SPM'; SPM(SPM <=0) = NaN;
    
    %latitude/long from the image
    lat = ncread(imageFile, 'lat'); lat = lat';
    lon = ncread(imageFile, 'lon'); lon = lon';
end

subaxis(5,4,[15 19],'Spacing',0.01,'MR',0.1);
m_proj('mercator','lon',  [-52.3 -51.8], 'lat',[-32.4 -31.7]);
m_pcolor(lon,lat,0.0993.*SPM.^(1.354));
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linewidth', 0.6, 'tickdir', 'none','gridcolor','none', 'backgroundcolor',[0.901960784313726 0.901960784313726 0.901960784313726],'fontsize',12);
m_text(-52.25,-32.35,'L8-OLI (2014-01-21)','color','w','fontsize',12);

caxis([0 100])

hold on

list_scenes = dir([pathDir 'L9*L2W.nc']);

for i=1
    
    fileName = char({list_scenes(i).name});
    imageFile = [pathDir fileName];
    SPM = ncread(imageFile,'TUR_Nechad2009_865');  SPM = SPM'; SPM(SPM <=0) = NaN;
    
    %latitude/long from the image
    lat = ncread(imageFile, 'lat'); lat = lat';
    lon = ncread(imageFile, 'lon'); lon = lon';
end

subaxis(5,4,[16 20],'Spacing',0.01,'MR',0.1);
m_proj('mercator','lon',  [-52.3 -51.8], 'lat',[-32.4 -31.7]);
m_pcolor(lon,lat,0.0993.*SPM.^(1.354));
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linewidth', 0.6, 'tickdir', 'none','gridcolor','none', 'backgroundcolor',[0.901960784313726 0.901960784313726 0.901960784313726],'fontsize',12);
m_text(-52.25,-32.35,'L9-OLI (2022-06-12)','color','k','fontsize',12);

caxis([0 100])

hold on

a = m_line(-52.10556,-32.02349,'marker','o','color',[0 0 0],'linewi',1,...
    'linest','none','markersize',8,'markerfacecolor',[0 0 0]);

c = m_line(-52.0954333,-32.2454167,'marker','<','color',[0 0 0],'linewi',1,...
    'linest','none','markersize',9,'markerfacecolor',[1 1 1]);

legend([a c],{'Buoy RS1', 'Buoy RS4'})
