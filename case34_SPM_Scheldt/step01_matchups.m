%------------------------------------------------------------------------ %
%                                                                         %
% MATLAB CODE TO READ SATELLITE SCENES (Sentinel 2A,B), FIND MATCHUPS     %
% WITH FIELD SPM (RWS) FOR WESTERSCHELDE ESTUARY                          %                                                     %
%                                                                         %
% THIS CODE REQUIRES:                                                     %
%                                                                         %
% Inputs:                                                                 %
% 1)Satellite L2 (atmospherically corrected) Sentinel 2                   %
% 2)folder SPM_RWS with in-situ SPM data + funcions to import data        %
%   ** make sure that the downloaded in-situ data is in UTC or correct time
%                 %
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

cd 'scheldt/SPM_RWS/'; %set folder path to data
files = dir('2022*.txt'); 

fileName = char({files(1).name});
[SPMfield,~,~,~,~] = process_fieldSPM(fileName);

weschelde_station = string(unique(SPMfield.Stations));

ll=1;
for i=1:length(weschelde_station)
    
    [st,~]=find(SPMfield.Stations==weschelde_station{(i)});
    SPMdata_st = SPMfield(st,:);
    SPMdata1 = removevars(SPMdata_st,'Stations');
    
    SPM_daily{:,ll} = SPMdata1;
    Lat_daily{:,ll} = SPMdata1.Lat;
    Lon_daily{:,ll} = SPMdata1.Lon;
    
    ll=ll+1;
end

%-------------------------------------------------------------------------%
%                      process satellite SPM data                         %
%-------------------------------------------------------------------------%

% directory structure pointing to where images are located
pathDir = '/Volumes/Expansion/Satellite_imagery/Scheldt/ACOLITE/Sentinel2/level2/';
list_scenes = dir([pathDir 'S2*L2W.nc']);


match_ups_S2 = [];

for i=1:length(list_scenes)
    
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
    
    satdate = datetime(str2double(fileName(9:12)),str2double(fileName(14:15)),str2double(fileName(17:18)),'Format','dd.MM.yyyy');
    
    sat_time = ncreadatt(imageFile,'/','isodate');
    sat_time = datetime(sat_time(12:19), 'Format', 'HH:mm:ss');
    
    
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
%                  find match-up data and collect data                    %
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
    
    
    for nn = 1:length(weschelde_station)
        
        SPM_field = SPM_daily{:,nn};
        
        [x,l] = find( (year(satdate(:,1))  == year(SPM_field.datetime)) & ...
                      (month(satdate(:,1)) == month(SPM_field.datetime)) & ...
                      (day(satdate(:,1))   == day(SPM_field.datetime)));
        
        
        field_time      = SPM_field.time(x); field_time.Format = 'HH:mm:ss';
        field_time.Hour = field_time.Hour - 1; %correct time difference (CET/MET to UTC)
        
        time_dif = abs((sat_time.Hour.*60 + sat_time.Minute) - ...
                     (field_time.Hour.*60 + field_time.Minute));
        
        
        if any(time_dif <= 30) %time interval 
            dataSPM = SPM_field.SPM(x,:);      dataSPM(time_dif>30)=[];
            LAT     = SPM_field.Lat(x,:);      LAT(time_dif>30)=[];
            LON     = SPM_field.Lon(x,:);      LON(time_dif>30)=[];
            dateSPM = SPM_field.datetime(x,:); dateSPM(time_dif>30)=[];
                                               field_time(time_dif>30)=[];
        else LAT = NaN;
        end
        
        if any(~isnan(LAT)) && any(x>0)
            
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
            

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
%                        create table with data                           %
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
            
            % TABELAO
            tabledata    = array2table([round(SPM_sample,2),round(SPMstd,2),lat_sample',lon_sample'],...
                              'VariableNames',{'SPM_nechad','SPM_nechad_std','lat_sample','lon_sample'});
            fielddata    = array2table([round(LAT,3), round(LON,3), ...
                              round(dataSPM,2)],'VariableNames',...
                              {'fieldLAT','fieldLON','field_SPM'});
            date         = table(datetime(datestr(dateSPM),'Format','dd.MM.yyyy'),...
                              repmat(datetime(sat_time,'Format','HH:mm:ss'),size(dateSPM,1),1),...
                              datetime(field_time,'Format','HH:mm:ss'),...
                              'VariableNames',{'dd.mm.yyyy','sat_HH.mm.ss_UTC',...
                              'field_HH.mm.ss_UTC'});
            st           = table(weschelde_station(nn),'VariableNames',{'Station'});
            
            match_ups_S2 = [match_ups_S2; date, st, fielddata, tabledata];
            
        end
    end
    
    disp(i)
    
    clear x
end

match_ups_S2 = sortrows(match_ups_S2,1,'ascend'); %sorting sampling 
match_ups_S2(isnan(match_ups_S2.SPM_nechad),:) = []; % removing NAN retrievals


%-------------------------------------------------------------------------%
%                           plotting match-ups                            %
%-------------------------------------------------------------------------%


%determine in-situ stations in matchups to plot
stations = unique(match_ups_S2.Station);
symb={'o','>', '<', 's'};

figure(1)
for i=1:length(stations)
    idx = ismember(match_ups_S2.Station,stations(i));
    hold all
    errorbar( match_ups_S2.field_SPM(idx), match_ups_S2.SPM_nechad(idx), ...
       match_ups_S2.SPM_nechad_std(idx) ,'k', 'LineStyle','none'); 
    
    dataplot(i) = scatter(match_ups_S2.field_SPM(idx), match_ups_S2.SPM_nechad(idx),20,...
                  symb{i},'MarkerFaceColor', 'k' ,'MarkerEdgeColor', 'k' ,'LineWidth',2);
end

oneone = plot(([0.01:100]),([0.01:100]),':k');
set(gca,'FontSize',12)
xlim([0 100])
ylim([0 100])
xlabel('SPM in-situ [gm^{-3}]','FontSize',12)
ylabel('SPM derived [gm^{-3}]','FontSize',12)
set(gca,'color',[1 1 1]); %sets background color
box on 
