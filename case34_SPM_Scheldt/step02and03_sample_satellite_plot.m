%------------------------------------------------------------------------ %
%                                                                         %
% MATLAB CODE TO READ SATELLITE SCENES (Sentinel 2A,B), sample SPM and    %
% compare with in-situ data                                               %  
%                                                                         %
%       Here we plot information relative to step 02 and step 03          %
%                                                                         %
% THIS CODE REQUIRES:                                                     %
%                                                                         %
% Inputs:                                                                 %
% 1)Satellite L2 (atmospherically corrected) Sentinel 2                   %
% 2)folder SPM_RWS with in-situ SPM data + funcions to import data        %
%                                                                         %
%                                                                         %
%                                                                         %
% developed by:                                                           %
% Juliana Tavora (j.tavora@utwente.nl)                                    %
% University of Twente                                                    %
% version November 2022                                                   %
%                                                                         %
%------------------------------------------------------------------------ %

close all
clear all
clc

%-------------------------------------------------------------------------%
%                        process in-situ SPM data                         %
%-------------------------------------------------------------------------%

files = dir('.txt');

fileName = char({files(1).name});
[SPMfield,~,~,~,~] = process_fieldSPM(fileName);

%%    ---------------------------   SPM    ----------------------------   %%

% directory structure pointing to where images are located
pathDir = 'pathto/ACOLITE/Sentinel2/level2/';
list_scenes = dir([pathDir 'S2*L2W.nc']);


timeseries_S2 = [];

LAT = SPMfield.Lat(1);
LON = SPMfield.Lon(1);


for i=1:length(list_scenes)
    
    fileName = char({list_scenes(i).name});
    imageFile = [pathDir fileName];
    
    satdate = datetime(str2double(fileName(9:12)),str2double(fileName(14:15)),str2double(fileName(17:18)),'Format','dd.MM.yyyy');
        
    
    %latitude\long from the image
    lat = ncread(imageFile, 'lat'); lat = lat';
    lon = ncread(imageFile, 'lon'); lon = lon';
    
    if fileName(3)=='A'
        SPM = ncread(imageFile,'SPM_Nechad2016_865');  SPM = SPM'; SPM(SPM ==0) = NaN;
    else
        SPM = ncread(imageFile,'SPM_Nechad2016_864');  SPM = SPM'; SPM(SPM ==0) = NaN;
    end
    
    
    minDist = [abs(lon - LON) + abs(lat - LAT)];
    [pLat, pLon] = find(minDist == min(abs(minDist(:))));
    
    
    %x size as array of begin\end points for box
    lon1 = pLon - 1;lon2 = pLon + 1;
    lat1 = pLat - 1;lat2 = pLat + 1;
    
    
    if  ~isnan(pLat)
        % find lat/lon on the scene
        minDist = abs(lon - LON) + abs(lat - LAT);
        [pLat, pLon] = find(minDist == min(abs(minDist(:))));
        
        % set sample box size (box size here = 3) as array of begin/end points for box
        lon1 = pLon - 1; lon2 = pLon + 1;
        lat1 = pLat - 1; lat2 = pLat + 1;
        
        if  ~isnan(pLat)
            SPM_sample      = nanmedian(SPM(lat1:lat2,lon1:lon2),'all');
            SPMstd          = std(SPM(lat1:lat2,lon1:lon2),0,'all','omitnan');
            lat_sample      = nanmedian(lat(lat1:lat2,lon1:lon2),'all');
            lon_sample      = nanmedian(lon(lat1:lat2,lon1:lon2),'all');
        end
        
        % TABELAO
        
        tabledata     = array2table([round(SPM_sample',2), round(SPMstd',2),...
            round(lon_sample',3), round(lat_sample',3)],...
            'VariableNames',{'SPM','SPM_std','lon_sample','lat_sample'});
        date          = datestr(satdate);
        date          = table(datetime(date,'Format','dd.MM.yyyy'),'VariableNames',{'dd.mm.yyyy'});
        sat_scene     = array2table(string(fileName), 'VariableNames',{'Scene'});
        timeseries_S2 = [timeseries_S2; date, sat_scene,tabledata];
    end
    
    disp(i)
    
end

timeseries_S2(isnan(timeseries_S2.SPM),:) =[];
%save('Sentinel2_timeseries.mat','timeseries_S2')


%-------------------------------------------------------------------------%
%                        plot figure of step 02                           %
%-------------------------------------------------------------------------%

figure(1);

%plot in-situ data
data = SPMfield.SPM;
field = histogram(log10(data),25,'Normalization','probability','FaceColor','k');
edge = field.BinEdges;
hold on
set(gca,'FontSize',12)

data = timeseries_S2.SPM;
sat = histogram(log10(data),edge,'Normalization','probability','FaceColor',[0.47,0.67,0.19]);
hold on
set(gca,'FontSize',12)


%-------------------------------------------------------------------------%
%                        plot figure of step 03                           %
%-------------------------------------------------------------------------%
figure(2);

%plot in-situ data
plot(SPMfield.datetime,SPMfield.SPM,'o')
xlim([datetime("2015-01-01") datetime("2020-12-31")])
hold on
plot(timeseries_S2.("dd.mm.yyyy"),timeseries_S2.SPM,'>')

