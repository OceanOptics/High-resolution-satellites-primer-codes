%-------------------------------------------------------------------------%
%                 plotting study case Coast of Maine - USA                %
%                                                                         %
%  inputs:                                                                %
%       - Sentinel 2A scene:  '20210618T153716_20mpolymer_S2A_OCX_chla_spm.nc'
%       - in-situ and satelite sampled data at in-situ location:   'dmc_matchup20162020.mat'
%       - m_map package                                                   %
%       - MBARI type regression code                                      %
%                                                                         %
%  output:                                                                %
%       - figure 6 of primer                                              %
%                                                                         %
%                                                                         %
% by: Juliana Tavora;  Binin Jiang.                                       %
% email: j.tavora@utwente.nl; lanzhiqishi@gmail.com                       %
%                                                                         %
% November 2022                                                           %
%                                                                         %
%-------------------------------------------------------------------------%


clear all
close all
clc


% ---------------------          read nc file         --------------------%
pathDir = 'path_to/';
list_scenes = dir([pathDir '*.nc']);
fileName = char({list_scenes(1).name});
imageFile = [pathDir fileName];
Chla = ncread(imageFile,'Chla');  
Chla = rot90(Chla') ;

% set lat/lon limits for satellite scene
long = linspace(-70-15/60-10/3600,-68-52/60-42/3600,5490);
lati = linspace(43+15/60+30/3600,  44+15/60+12/3600,5490);
[lon,lat] = meshgrid(long,lati);

% plot panel d
figure(1)
subaxis(2,3,[3 6],'Spacing',0.01,'MR',0.1);
m_proj('mercator','lon',  [-69-39/60-00/3600 -69-27/60-00/3600], 'lat',[43+40/60+00/3600  44+4/60+00/3600]);
m_pcolor(lon+0.002,lat,medfilt2(Chla,[2, 2]));
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linewidth', 0.6, 'tickdir', 'none','gridcolor','none', 'backgroundcolor',[0.901960784313726 0.901960784313726 0.901960784313726],'fontsize',12);

a = m_line(-69.579562,43.934180,'marker','o','color',[0 0 0],'linewi',1,...
    'linest','none','markersize',8,'markerfacecolor',[0 0 0]);

caxis([0 5])


% load in-aitu and satelite data
load('dmc_matchup20162020.mat')

%plot panel c
subaxis(2,3,[4 5],'Spacing',0.01,'MR',0.1);

%plot background
for i=[2016:1:2022]; yy(i-2015) = datetime([i 1 15]); end
for i=2:2:6
    fill([yy(i)  yy(i) ...
        yy(i+1) yy(i+1)],...
        [0 599 599 0], [0.98 0.98 0.98],'EdgeColor','none');
    hold on
end
clear yy i 

hold on
a = scatter(DMC_time_20162020,DMC_chla_20162020,'o','MarkerFaceColor','w','MarkerEdgeColor','k');
b = plot(sentineltime_realdate_order,OCX_sh_median_all,'>','MarkerFaceColor','k','MarkerEdgeColor','k');
xlim([datetime("2016-01-01") datetime("2020-12-31")])
ylim([0 15])
set(gca,'FontSize',12)

ylabel('Sentinel2  chla  ($\mu$g/l)', 'interpreter', 'tex')
xlabel('Year')
legend([a b],{'In-situ Chla','Sentinel 2A,B derived Chla'})


% find matchus between in-situ and satelite
s = length(DMC_time_20162020);
[m,n]=size(sentineltime_realdate_order);

for i=1:m
    
    time_DMC_chla      =     find(abs(datenum(DMC_time_20162020) - datenum(sentineltime_realdate_order(i,1))) < 1e-1);%find the location
    
    DMC_Chla_sentinel2 =     DMC_chla_20162020(time_DMC_chla);
    DMC_Chla_sentinel2(isempty(DMC_Chla_sentinel2)== 1 ) =   NaN;
    
    eval(['chla_DMCBuoy_mean_chlaall(i)','=','DMC_Chla_sentinel2(1,1)',';']); 
end


subaxis(2,3,[1],'Spacing',0.01,'MR',0.1);
hold all
errorbar(chla_DMCBuoy_mean_chlaall,OCX_sh_median_all,OCX_sh_std_all,'Horizontal','k', 'LineStyle','none'); %percentual
dataplot = scatter(chla_DMCBuoy_mean_chlaall,OCX_sh_median_all,'MarkerFaceColor', 'w' ,'MarkerEdgeColor', 'k' ,'LineWidth',1);


%type 2 regression
[xData, yData] = prepareCurveData(chla_DMCBuoy_mean_chlaall,OCX_sh_median_all);
[m,b,~,~,~]  = lsqbisec(xData,yData);
latlot       = linspace(0.001,400,500);
lonlot       = m.*latlot + b;

hold on
oneone = plot(([0.001:400]),([0.001:400]),':k');
reg= plot(latlot,lonlot,'k','LineWidth',0.5);
set(gca,'FontSize',12)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([0.1 50])
ylim([0.1 50])
xlabel('Sentinel 2 derived Chla  ($\mu$g/l)', 'interpreter', 'latex','FontSize',12)
ylabel('In-situ Chla  ($\mu$g/l)', 'interpreter', 'latex','FontSize',12)
set(gca,'color',[1 1 1]); %sets grey background

legend([oneone reg],{'1:1','y = -0.97x + 3.97, n = 20'})

%metrics
chla_err        =  chla_DMCBuoy_mean_chlaall - OCX_sh_median_all;
chla_sam_number = numel(find(~isnan(chla_err)));
chla_error      = real(sqrt(nansum((chla_err).^2)./chla_sam_number));
RMSE_chla       = chla_error;
log_err         = log10( chla_DMCBuoy_mean_chlaall)- log10(OCX_sh_median_all);
MAPE_chla       = nanmedian(abs(chla_err)./ chla_DMCBuoy_mean_chlaall) .*100;
MAE_chla        = 10 .^(nansum(abs(log_err)) ./chla_sam_number );


%plot panel b
subaxis(2,3,[2],'Spacing',0.01,'MR',0.1);
idx = isnan(chla_DMCBuoy_mean_chlaall);
field = histogram(log10(chla_DMCBuoy_mean_chlaall(~idx)),25,'Normalization','probability','FaceColor','k');
hold on
set(gca,'FontSize',12)
ylim([0 0.3])

edge = field.BinEdges;
idx = isnan(OCX_sh_median_all);
sat = histogram(log10(OCX_sh_median_all(~idx)),edge,'Normalization','probability','FaceColor','w');
hold on
set(gca,'FontSize',12)


leg = legend([field, sat],{'in-situ: 28.1 $\pm$ 44.45,','satellite-derived: 6.26 $\pm$ 44.39'},'EdgeColor','none','interpreter','latex');
title(leg,' Turbidity [FTU]')

xlabel('Chla  [log_{10}($\mu$g/l)]', 'interpreter', 'latex','FontSize',12)
ylabel('Relative frequency','FontSize',12)





 
 
 


