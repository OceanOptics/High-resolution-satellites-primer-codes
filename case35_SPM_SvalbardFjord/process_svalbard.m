clear all
clc

load('matchups_acolite_collection2_level2.mat')
load('matchups_usgs_collection2_level2.mat')


%--------------------------- plot match-up  ------------------------------%
id = isnan(match_ups.SPM);
match_ups(isnan(match_ups.SPM),:)=[];

time_dif = ((match_ups.("sat_HH.mm.ss_UTC").Hour.*60 + match_ups.("sat_HH.mm.ss_UTC").Minute) - ...
    (match_ups.time.Hour.*60 + match_ups.time.Minute))./60; time_dif = time_dif(18:end);

figure(1)
subaxis(3,2,1,'Spacing',0.01,'MR',0.1);
hold all
err         = errorbar(match_ups.SPM(18:end), match_ups.SPM_mw(18:end),        match_ups.SPM_std(18:end),  'k', 'LineStyle','none'); %percentual(1:43)
dataplot2   = scatter(match_ups.SPM(18:end), match_ups.SPM_mw(18:end),70,  'MarkerFaceColor', 'w','MarkerEdgeColor', 'k','Marker','o');
dataplot3   = scatter(match_ups_usgs_l2c2.SPM, match_ups_usgs_l2c2.SPM_mw,70,  'MarkerFaceColor', 'k','MarkerEdgeColor', 'k','Marker','o');

oneone      = plot(([0.01:400]),([0.01:400]),':k');

%type 2 regression
[xData, yData] = prepareCurveData(match_ups.SPM(18:end), match_ups.SPM_mw(18:end));
[m,b,~,~,~]  = lsqbisec(xData, yData);
latlot       = linspace(1,400);
lonlot       = m.*latlot + b;

%metrics
spm_err        =  xData - yData;
spm_sam_number = numel(find(~isnan(spm_err)));
spm_error      = real(sqrt(nansum((spm_err).^2)./spm_sam_number));
RMSE_spm       = spm_error;
log_err        = log10(xData)- log10(yData);
MAPE_spm       = nanmedian(abs(spm_err)./ xData) .*100;
MAE_spm        = 10 .^(nansum(abs(log_err)) ./spm_sam_number );



[r_kendall,p]   = corr(xData,yData,'Type','Kendall','Rows','pairwise')

reg= plot(latlot,lonlot,'k','LineWidth',0.5);
set(gca,'FontSize',12)
xticklabels({'10^{0}', '10^{1}', '10^{2}'})
yticklabels({'10^{0}', '10^{1}', '10^{2}'})
xlabel('SPM in-situ [gm^{-3}]','FontSize',12)
ylabel('SPM derived [gm^{-3}]','FontSize',12)
set(gca,'color',[1 1 1]); %sets grey background
box on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([10 400])
ylim([10 400])

legend([oneone reg dataplot2 dataplot3],{'1:1','y = 0.82x - 10.46, r = 0.44 (NS); n = 9',...
    '2015/08/14 - ACOLITE L2R', '2015/08/14 - USGS level 2'})


% M MAP
addpath ~/Documents/MATLAB/m_map2
rehash toolboxcache

load('maps_SPM_MWalg_scenes3_usgs.mat')
load('maps_SPM_MWalg_secene3_15aug.mat')


figure(1)
subaxis(3,2,3, 'Spacing', 0.05,'SpacingHoriz',0.05,'Padding', 0, 'Margin', 0.1);
m_proj('mercator','lon', [17.06 17.46] , 'lat', [78.38 78.47]);
m_pcolor(lon_lala-0.00016,lat_lala-0.0037,log10(SPM_alg)); shading flat;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linewidth', 0.6, 'tickdir', 'none','gridcolor','none', 'backgroundcolor',[0.901960784313726 0.901960784313726 0.901960784313726],'fontsize',12);
colorbar 
caxis([0 2.5])

subaxis(3,2,4, 'Spacing', 0.05,'SpacingHoriz',0.05,'Padding', 0, 'Margin', 0.1);
m_proj('mercator','lon', [17.06 17.46] , 'lat', [78.38 78.47]);
m_pcolor(lon_lala-0.00016,lat_lala-0.0037,SPM_alg_err); shading flat;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linewidth', 0.6, 'tickdir', 'none','gridcolor','none', 'backgroundcolor',[0.901960784313726 0.901960784313726 0.901960784313726],'fontsize',12);
caxis([0 100])
colorbar


subaxis(3,2,5, 'Spacing', 0.05,'SpacingHoriz',0.05,'Padding', 0, 'Margin', 0.1);
m_proj('mercator','lon', [17.06 17.46] , 'lat', [78.38 78.47]);
m_pcolor(lon_lala-0.00016,lat_lala-0.0037,log10(SPM_alg_usgs)); shading flat;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linewidth', 0.6, 'tickdir', 'none','gridcolor','none', 'backgroundcolor',[0.901960784313726 0.901960784313726 0.901960784313726],'fontsize',12);
colorbar
caxis([0 2.5])


subaxis(3,2,6, 'Spacing', 0.05,'SpacingHoriz',0.05,'Padding', 0, 'Margin', 0.1);
m_proj('mercator','lon', [17.06 17.46] , 'lat', [78.38 78.47]);
m_pcolor(lon_lala-0.00016,lat_lala-0.0037,SPM_alg_err_usgs); shading flat;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linewidth', 0.6, 'tickdir', 'none','gridcolor','none', 'backgroundcolor',[0.901960784313726 0.901960784313726 0.901960784313726],'fontsize',12);
caxis([0 100])
colorbar




