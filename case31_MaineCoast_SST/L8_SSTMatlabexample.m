
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract temperature from landsat level 2 products %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear variables

% Requirements:

% (1) A landsat level 2 product bundel downloaded from https://earthexplorer.usgs.gov/
% (2) Download the getLandsatL2_SST function from the Ocean Optics github...
% at https://github.com/OceanOptics/getLandsatL2_SST
% (3) Make the folling changes to Matlabs geotiffread function outlined in...
% https://cosmojiang.wordpress.com/2018/04/02/matlab-geotiffread-for-multiple-layers/

%% Extract the latitude, longitude, and temperature data

% set this path to where you have saved the Level 2 file
% The command 'pwd' in the commond window will list the current path if...
% navigate to it in the current folder window on the left

% Run the function
pathDir = 'path_to/GEE_SST/LC08_L2SP_011030_20220909_20220914_02_T1/';

    % Run the function
    [dt, lat, lon, temperature, fmask] = getLandsatL2_SST(pathDir, false, [0 99]);

% dt returns the date of the image
% lat returns a matrix of the latitudes
% lon returns a matrix of the longitudes
% temperature is the matrix of SSTs
% fmask is contines identifiers for masked out pixels (land, ice, cloud, ect.)

%% Plot the image

figure(1)
subaxis(3,3,[3 6],'Spacing',0.01,'MR',0.1);
m_proj('mercator','lon',  [-69-39/60-00/3600 -69-27/60-00/3600], 'lat',[43+40/60+00/3600  44+4/60+00/3600]);
m_pcolor(lon, lat, temperature); 
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linewidth', 0.6, 'tickdir', 'none','gridcolor','none', 'backgroundcolor',[0.901960784313726 0.901960784313726 0.901960784313726]);

caxis([17 22])
cb = colorbar('Location','south','FontSize',12);
cb.Label.String = 'T [^{o}C]';

set(gca,'FontSize',12)

%% Find the closest pixel to site of interest

% Define coordinates for site of interest
latRef = 43.99818; % latitude
lonRef = -69.54253; % Longitude

% Use the findClosestPixel function from the bottom of the code
[row, column, closestPixel] = findClosestPixel(lon, lat, lonRef, latRef);
% returns the row, column, and index of the closest pixel to the ...
% reference site in the lon and lat matrix from the satellite image

%% Plot to make sure the extracted data is the location you expect it to be

% Use the row and column of the cosest pixel to the site to subset the... 
% matrices to 200 pixels on either side
latCrop = lat(row-200:row+200, column-200:column+200);
lonCrop = lon(row-200:row+200, column-200:column+200);
tempCrop = temperature(row-200:row+200, column-200:column+200);

% Full image with site
figure(2)
sstPlot = pcolor(lon, lat, temperature);
set(sstPlot, 'EdgeColor', 'none');
colorbar;
colormap('jet');
hold on
plot(lonRef, latRef, 'o',...
    'LineWidth', 10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','m',...
    'MarkerSize', 8);

% Cropped image with site
figure(3)
sstPlot = pcolor(lonCrop, latCrop, tempCrop);
set(sstPlot, 'EdgeColor', 'none');
colorbar;
colormap('jet');
hold on
plot(lonRef, latRef, 'o',...
    'LineWidth', 10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','m',...
    'MarkerSize', 8);

% Extract window of data from point of interest

the neighbors function returns all the temperature data from the pixel of ...
interest and 5x5 window of pixels around it using the row and column...
information from the findClosestPixel function

allVal = neighbors(row, column, temperature);

Average temperature at the site
meanTemp = mean(allVal);
standard deviation of temperature at the site
stdTemp = std(allVal);

%% Exract all touching values

function allVal = neighbors(row, column, dataMatrix)

% data from find closest pixel
x = column;
y = row;

% 3x3 touching values
pix1  = dataMatrix(y, x);
pix2  = dataMatrix(y+1, x);
pix3  = dataMatrix(y+1, x+1);
pix4  = dataMatrix(y, x+1);
pix5  = dataMatrix(y+1, x-1);
pix6  = dataMatrix(y-1, x);
pix7  = dataMatrix(y-1, x-1);
pix8  = dataMatrix(y, x-1);
pix9  = dataMatrix(y+1, x-1);
% 5x5 matrix
pix10 = dataMatrix(y+2, x-2);
pix11 = dataMatrix(y+2, x-1);
pix12 = dataMatrix(y+2, x);
pix13 = dataMatrix(y+2, x+1);
pix14 = dataMatrix(y+2, x+2);
pix15 = dataMatrix(y+1, x+2);
pix16 = dataMatrix(y, x+2);
pix17 = dataMatrix(y-1, x+2);
pix18 = dataMatrix(y-2, x+2);
pix19 = dataMatrix(y-2, x+1);
pix20 = dataMatrix(y-2, x);
pix21 = dataMatrix(y-2, x-1);
pix22 = dataMatrix(y-2, x-2);
pix23 = dataMatrix(y-1, x-2);
pix24 = dataMatrix(y, x-2);
pix25 = dataMatrix(y+1, x-2);

allVal = [pix1, pix2, pix3, pix4, pix5, pix6, pix7, pix8, pix9, pix10,...
          pix11, pix12, pix13,pix14, pix15, pix16, pix17, pix18, pix19,...
          pix20, pix21, pix22, pix23, pix24, pix25]';
end


%% Function for finding the closest points

function [row, column, closest_pixel] = findClosestPixel(lon, lat, lonref, latref)

% Subtract refernce points from grids
minDist = abs(lon - lonref) + abs(lat - latref);

% Find the index of the closest point
closest_pixel = find(minDist == min(abs(minDist(:))));

% return the row and column
[row, column] = find(minDist == min(abs(minDist(:))));

% site = [x, y];
% 
% row = site(1);
% column = site(2);

end

%%%%%%%%%%%%%%%%%%% Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code to reproduce Figure 3 

% Load data
load('landsatDREcomparison.mat')
load('landsat_GEE_all.mat')
load('landsat_GEE_90PercentClear.mat')
load('histogramData.mat')



subaxis(3,3,[1 2],'Spacing',0.01,'MR',0.1);


hold on
id = find(landsatGEEall.Sensor == 5);
doy = day(landsatGEEall.DateTime(id),'dayofyear');

% Panel a
a = scatter(doy,landsatGEEall.SST_C(id),'o','MarkerFaceColor','w','MarkerEdgeColor','k');

id = find(landsatGEEall.Sensor == 8);
doy = day(landsatGEEall.DateTime(id),'dayofyear');

% Panel c
c = scatter(doy,landsatGEEall.SST_C(id),'square','MarkerFaceColor',[0.9 0.9 0.9],'MarkerEdgeColor','k');

id = find(landsatGEEall.Sensor == 7);
doy = day(landsatGEEall.DateTime(id),'dayofyear');

% Panel b
b = scatter(doy,landsatGEEall.SST_C(id),'>','MarkerFaceColor','k','MarkerEdgeColor','k');

id = find(landsatGEEall.Sensor == 9);
doy = day(landsatGEEall.DateTime(id),'dayofyear');

% Panel d
d = scatter(doy,landsatGEEall.SST_C(id),'d','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k');


xlim([0 366])
ylim([-10 30])

set(gca,'FontSize',12)

ylabel('Sea Surface Temperature  (^oC)', 'interpreter', 'tex')
xlabel('Day of Year')
legend([a b c d],{'L5','L7','L8','L9'})



subaxis(3,3,[4 5],'Spacing',0.01,'MR',0.1);


hold on
id = find(landsatGEE90PercentClear.Sensor == 5);
doy = day(landsatGEE90PercentClear.DateTime(id),'dayofyear');

a = scatter(doy,landsatGEE90PercentClear.SST_C(id),'o','MarkerFaceColor','w','MarkerEdgeColor','k');

id = find(landsatGEE90PercentClear.Sensor == 8);
doy = day(landsatGEE90PercentClear.DateTime(id),'dayofyear');

c = scatter(doy,landsatGEE90PercentClear.SST_C(id),'square','MarkerFaceColor',[0.9 0.9 0.9],'MarkerEdgeColor','k');

id = find(landsatGEE90PercentClear.Sensor == 7);
doy = day(landsatGEE90PercentClear.DateTime(id),'dayofyear');

b = scatter(doy,landsatGEE90PercentClear.SST_C(id),'>','MarkerFaceColor','k','MarkerEdgeColor','k');

id = find(landsatGEE90PercentClear.Sensor == 9);
doy = day(landsatGEE90PercentClear.DateTime(id),'dayofyear');

d = scatter(doy,landsatGEE90PercentClear.SST_C(id),'d','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k');


xlim([0 366])
ylim([-10 30])

set(gca,'FontSize',12)

ylabel('Sea Surface Temperature  (^oC)', 'interpreter', 'tex')
xlabel('Day of Year')

subaxis(3,3,7,'Spacing',0.01,'MR',0.1);

a = scatter(landsatDREcomparison.Lobo_sst,landsatDREcomparison.Landsat_sst,'o','MarkerFaceColor','w','MarkerEdgeColor','k');
hold on
oneone = plot(([2:0.1:29]),([2:0.1:29]),':k');
xlim([2 29]); ylim([2 29])

%type 2 regression
[xData, yData] = prepareCurveData(landsatDREcomparison.Lobo_sst,landsatDREcomparison.Landsat_sst);
[m,b,~,~,~]  = lsqbisec(xData,yData);
latlot       = linspace(0.001,400,500);
lonlot       = m.*latlot + b;

reg= plot(latlot,lonlot,'k','LineWidth',0.5);

box on

ylabel('SST GEE  (^oC)', 'interpreter', 'tex')
xlabel('SST in-situ (^oC)')


subaxis(3,3,8,'Spacing',0.01,'MR',0.1);
situ = histogramData(57:end,:);
sat = histogramData(1:56,:);

[ida,idb] = ismember(sat.Date,situ.Date);
idb(idb<1) =[]

a = scatter(situ.Temp_C(idb),sat.Temp_C(ida),'o','MarkerFaceColor','w','MarkerEdgeColor','k');
hold on
oneone = plot(([2:0.1:29]),([2:0.1:29]),':k');
xlim([2 29]); ylim([2 29])

%type 2 regression
[xData, yData] = prepareCurveData(situ.Temp_C(idb),sat.Temp_C(ida));
[m,b,~,~,~]  = lsqbisec(xData,yData);
latlot       = linspace(0.001,400,500);
lonlot       = m.*latlot + b;

reg= plot(latlot,lonlot,'k','LineWidth',0.5);

box on

ylabel('SST downloaded Landsat  (^oC)', 'interpreter', 'tex')
xlabel('SST in-situ (^oC)')


subaxis(3,3,9,'Spacing',0.01,'MR',0.1);

field = histogram((histogramData.Temp_C(57:end)),20,'Normalization','probability','FaceColor','k');
hold on
set(gca,'FontSize',12)
ylim([0 0.2])

edge = field.BinEdges;
sat = histogram((histogramData.Temp_C(1:56)),edge,'Normalization','probability','FaceColor','w');
hold on
set(gca,'FontSize',12)


leg = legend([field, sat],{'in-situ: 15.67 \pm 5.09,','satellite-derived: 15.34 \pm 5.70'},'EdgeColor','none','interpreter','tex');

ylabel('Relative Frequency', 'interpreter', 'tex')
xlabel('SST (^oC)')
