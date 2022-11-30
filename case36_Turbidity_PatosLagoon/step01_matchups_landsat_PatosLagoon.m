%-------------------------------------------------------------------------- %
%                                                                           %
% MATLAB CODE TO READ SATELLITE SCENES (L5, L7, L8 and L9), FIND MATCHUPS   %
% WITH FIELD SPM (from SIMCOSTA) FOR Patos Lagoon ESTUARY                   %                                                 
%                                                                           %
% THIS CODE REQUIRES:                                                       %
%                                                                           %
% Inputs:                                                                   %
% 1)Satellite Landsat *L2W.nc file (atmospherically corrected using Acolite)%
% 2)folder in-situ Turbidity data + funcions to import data                 %
%   ** make sure that the downloaded in-situ data is in UTC or correct time %
%   ** you can use MATLABs 'Import Data' optin the Menu bar                 %
%                                                                           %
%                                                                           %
% Outputs:
% 1) table containing match-up turidity calulated from Nechad et al., 2009 algorithm 
%                                                                           %
% developed by:                                                             %
% Juliana Tavora (j.tavora@utwente.nl)                                      %
% University of Twente                                                      %
% version November 2022                                                     %
%                                                                           %
%-------------------------------------------------------------------------- %

clear all
close all
clc

cd 'X:\LP_Processadas\L8\'


%-------------------------------------------------------------------------%
%                      process satellite SPM data                         %
%-------------------------------------------------------------------------%



% directory structure pointing to where images are located
mainDir = '/Volumes/JTavora_ITC/PatosLagoon/';


%mainDir = pathDir;


list_scenes = dir([mainDir 'L*_L2W.nc']);  % get list of all .nc files
x = size(list_scenes);
numFiles = x(1); % get number of files found


big_table_matchup = [];

for i=28: size(list_scenes,1)
    
    x2 = {list_scenes(i).name};
    B = convertStringsToChars(string(x2)');
    
    satdate = datetime(str2double(B(8:11)),str2double(B(13:14)),str2double(B(16:17)), ...
        'Format','dd-MMM-yyyy')
    
    %read in ALL files in the directory, retrieving the data
    imageFile = [pathDir char(x2)];
    
    sat_time = ncreadatt(imageFile,'/','isodate');
    sat_time = datetime(sat_time(12:19), 'Format', 'HH:mm:ss');
    
    %latitude/long from the image
    lat = ncread(imageFile, 'lat');
    lon = ncread(imageFile, 'lon');
    
    %---------------------------------------------------------------------%
    % since earliest in-situ Turbidity data was collected in 2016, the only
    % satellite/sensor to be sampled is either L8 or L9, hence: 
    
    if str2num(B(2)) == 8
        %rhos are already corrected by acolite
        Turb_Nechad16_865 = ncread(imageFile, 'TUR_Nechad2016_865');  Turb_Nechad16_865(Turb_Nechad16_865 <0 ) = NaN;
    else
        %rhos are already corrected by acolite
        Turb_Nechad16_865 = ncread(imageFile, 'TUR_Nechad2016_865');  Turb_Nechad16_865(Turb_Nechad16_865 <0 ) = NaN;
    end
    

    Turb_Nechad16_865 = Turb_Nechad16_865';    
    lat =  lat';
    lon = lon';
 
    
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
    %                  find match-up data and collect data                    %
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
    LON = [-52.10556,-52.0980333,-52.0954333]; 
    LAT = [-32.02349,-32.1345833,-32.2454167];
    
    for ii=1:length(LAT)
        
        data = Turb_boia_data{:,ii};
        [x,l] = find( (year(satdate(:,1)) == year(data.datetime)) & ...
            (month(satdate(:,1)) == month(data.datetime)) & ...
            (day(satdate(:,1)) == day(data.datetime)));
        
        
        field_time = data.datetime(x); field_time.Format = 'HH:mm:ss';
        
        time_dif = abs((sat_time.Hour.*60 + sat_time.Minute) - ...
            (field_time.Hour.*60 + field_time.Minute));
        
        [~,r_min] = (nanmin(time_dif));
        
        time_dif=time_dif(r_min); x = x(r_min);
        
        clear nodata
        
        if any(time_dif <= 30)
            dataTurb = data(x,:);
            field_time=field_time(r_min);
            dateT = dataTurb.datetime;
        else
            nodata = 1; %LAT(ii) = NaN;
        end
        
        % find 'address' of lat/lon from in-situ data
        if exist('nodata','var') ==0 %&& any(x>0)
            
            % find lat/lon on the scene
            minDist = [abs(lon - LON(ii)) + abs(lat - LAT(ii))];
            [pLat, pLon] = find(minDist == min(abs(minDist(:))));
            
            % set sample box size (box size here = 5) as array of begin/end points for box
            lon1 = pLon - 2; lon2 = pLon + 2;
            lat1 = pLat - 2; lat2 = pLat + 2;
            
            if  ~isnan(pLat)
                Turb_Nechad16_865_box_median = squeeze(nanmedian(Turb_Nechad16_865(lat1:lat2,lon1:lon2),[1 2]))';
                Turb_Nechad16_865_box_std    = squeeze(std(Turb_Nechad16_865(lat1:lat2,lon1:lon2),0,'all','omitnan'))';
                
                lat_box_median = squeeze(nanmedian(lat(lat1:lat2,lon1:lon2),[1 2]))';
                lon_box_median = squeeze(nanmedian(lon(lat1:lat2,lon1:lon2),[1 2]))';
            end
            
            
            %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
            %                        create table with data                           %
            %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
            
            tabledata    = array2table([Turb_Nechad16_865_box_median, Turb_Nechad16_865_box_std, ...
                round(lon_box_median',3),round(lat_box_median',3)],...
                'VariableNames',{'Turb_Nechad16_865_median','Turb_Nechad16_865_std',...
                'lon_sample','lat_sample'});
            latlondata    = array2table([round(LAT(ii),3), round(LON(ii),3)],'VariableNames',...
                {'fieldLAT','fieldLON'});
            date          = table(datetime(datestr(dateT),'Format','dd.MM.yyyy'),...
                datetime(sat_time,'Format','HH:mm:ss'),...
                datetime(field_time,'Format','HH:mm:ss'),...
                'VariableNames',{'dd.mm.yyyy','sat_HH.mm.ss_UTC',...
                'field_HH.mm.ss_UTC'});
            
            big_table_matchup = [big_table_matchup; date, dataTurb,tabledata];
            
            clear nodata x 
        end
    end
end

big_table_matchup(isnan(big_table_matchup.Turb_Nechad16_865_median),:) = [];

save('Landsat_matchups.mat','big_table_matchup')
