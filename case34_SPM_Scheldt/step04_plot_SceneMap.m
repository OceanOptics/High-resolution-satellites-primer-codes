%------------------------------------------------------------------------ %
%                                                                         %
%     MATLAB example CODE TO READ SATELLITE SCENES (Sentinel 2A,B) and    %                                                     %
%                              plot scenes                                %
%
% THIS CODE REQUIRES:                                                     %
%                                                                         %
% Inputs:                                                                 %
% 1) Satellite L2 (atmospherically corrected) Sentinel 2                  %
% 2) m_map package added to toolbox of matlab                             %
%                                                                         %
%  - Here the user can adjust the code to plot all scenes, as wished.     %
%  - m-map package has more funtionalities than the ones presented here,  %
%    for more please visit the website                                    %
%                                                                         %
% developed by:                                                           %
% Juliana Tavora (j.tavora@utwente.nl)                                    %
% University of Twente                                                    %
% version November 2022                                                   %
%                                                                         %
%------------------------------------------------------------------------ %

%directory of scenes
pathDir = 'D:\S2_Scheldt\level2\';
list_scenes = dir([pathDir 'S2*L2W.nc']);


timeseries_S2 = [];
for i=1% :length(list_scenes)
    
    fileName = char({list_scenes(i).name});
    imageFile = [pathDir fileName];
       
    if fileName(3)=='A'
        SPM = ncread(imageFile,'SPM_Nechad2016_865');  SPM = SPM'; SPM(SPM ==0) = NaN;
    else
        SPM = ncread(imageFile,'SPM_Nechad2016_864');  SPM = SPM'; SPM(SPM ==0) = NaN;
    end
    
    %latitude/long from the image
    lat = ncread(imageFile, 'lat'); lat = lat';
    lon = ncread(imageFile, 'lon'); lon = lon';
   
end


figure(1)
m_proj('mercator','lon',  [3.19 4.42], 'lat',[51.15 51.58]);
m_pcolor(lon,lat,SPM); 
m_grid('linewidth', 0.6, 'tickdir', 'none','gridcolor','none', 'backgroundcolor',[0.901960784313726 0.901960784313726 0.901960784313726]);
m_grid('linestyle','none','box','fancy','tickdir','out');

caxis([0 250])
cb = colorbar('Location','south','FontSize',12);
cb.Label.String = 'SPM [g.m^{-3}]';

hold on

% optional to plot location of in-situ station
% declare LON and LAT
% m_line(LON,LAT,'marker','o','color',[0.929411764705882 0.694117647058824 0.125490196078431],'linewi',1,...
%    'linest','none','markersize',8,'markerfacecolor',[0.929411764705882 0.694117647058824 0.125490196078431]);
