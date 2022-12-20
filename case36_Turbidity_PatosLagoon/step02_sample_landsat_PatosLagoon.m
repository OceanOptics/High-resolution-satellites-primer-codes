
%------------------------------------------------------------------------ %
%                                                                         %
% MATLAB CODE TO PLOT data from SATELLITE SCENES (L5, L7, L8 and L9) from
% Patos Lagoon ESTUARY                 %                                                 
%                                                                         %
% THIS CODE REQUIRES:                                                     %
%                                                                         %
% Inputs:                                                                 %
% 1) folder with in-situ Turbidity data and import file function          %
% 2) L5 scenes                                                            %
% 3) L7 scenes                                                            %
% 4) L8 scenes                                                            %
% 5) L9 scenes                                                            %
%  
% Outputs:
% 1) .mat file with all Buoy RS4 data
% 2) .mat file with turbidity calculated from all L5, L7, L8, L9 scenes
%
%
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



%-------------------------------------------------------------------------%
%                     process in-situ Turibidity data                     %
%-------------------------------------------------------------------------%

cd 'path_to/Turbidity_SIMCOSTA/'%set directory where in-situ Turbidity and import function are

% getting RS4 (buoy 4 only) 
list_files = dir(['SIMCOSTA_RS-4*.csv']);  % get list of all .csv turbidity files in directory
x = size(list_files);


boias = {"RS4"};
Turb_boias = [];
for i =1:x
    % boia RS1
    Turb = import_T_SIMCOSTA(list_files(i).name); %RS1=30946 %RS2=43986 %RS4= 27283
    Turb_date = datetime([Turb.YEAR Turb.MONTH Turb.DAY Turb.HOUR Turb.MINUTE Turb.SECOND]);
    Turb_date = array2table(Turb_date);
    Turb(:,8) = Turb_date;
    Turb(:,1:6) =[];
    
    a = convertStringsToChars(boias{i}); 
    if a(3)==list_files(i).name(13)
        Turb(:,3) = array2table(repmat(boias{i},height(Turb),1));
    end
    Turb.Properties.VariableNames{2} = 'datetime'; clear Turb_date
    
    Turb_buoy4 = [Turb_buoy4; Turb]; clear Turb
end
Turb_buoy4.Properties.VariableNames{3} = 'boias';
Turb_buoy4(isnan(Turb_buoy4.Avg_Turb),:)=[];

save('Turb_Buoy4.mat','Turb_buoy4')

%-------------------------------------------------------------------------%
%                    process satellite turbidity data                     %
%-------------------------------------------------------------------------%

%  directory structure pointing to where images are located
mainDir = '/path_to/';
list_scenes = dir([mainDir 'L*_L2W.nc']);  % get list of all .nc files


% lat/ lon of insitu station relativo de Buoy RS4
LAT = -32.2454167;	
LON = -52.0954333;	

big_table_Landsat_RS4 = [];

for i=size(list_scenes,1)
    
    x2 = {list_scenes(i).name};
    B = convertStringsToChars(string(x2)');
    
    if B(2) == num2str(5)
        day     = str2double(B(15:16));
        month   = str2double(B(12:13));
        year    = str2double(B(7:10));
    else
        day     = str2double(B(16:17));
        month   = str2double(B(13:14));
        year    = str2double(B(8:11));
    end
    
    
    satdate = datetime(year,month,day, 'Format','dd.MM.yyyy');

    % read in ALL files in the directory, retrieving the data
    imageFile = [mainDir char(x2)];   
    
    %latitude\long from the image
    lat = ncread(imageFile, 'lat'); lat = lat';
    lon = ncread(imageFile, 'lon'); lon = lon';
    
    %NIR bands of landsat sensors differ
    if     B(2) == num2str(5)
                Turb = ncread(imageFile,'TUR_Nechad2009_839');  Turb = Turb'; Turb(Turb <=0) = NaN;
    elseif B(2) == num2str(7)
                Turb = ncread(imageFile,'TUR_Nechad2009_835');  Turb = Turb'; Turb(Turb <=0) = NaN;
    elseif B(2) == num2str(8) | B(2) == num2str(9)
                Turb = ncread(imageFile,'TUR_Nechad2009_865');  Turb = Turb'; Turb(Turb <=0) = NaN;      
    end
    
    flag =  ncread(imageFile, 'l2_flags');  flag = double(flag); flag(flag >0) = NaN; flag(flag <1) = 1;
    
    Turb = Turb.*flag';
    
    
    minDist = [abs(lon - LON) + abs(lat - LAT)];
    [pLat, pLon] = find(minDist == min(abs(minDist(:))));
    
    %x size as array of begin\end points for box
    lon1 = pLon - 1;lon2 = pLon + 1;
    lat1 = pLat - 1;lat2 = pLat + 1;
    
    
    if  ~isnan(pLat)
        Turb_NIRsample  = Turb(lat1:lat2,lon1:lon2);   % get the sample box
        latsample       = lat(lat1:lat2,lon1:lon2);
        lonsample       = lon(lat1:lat2,lon1:lon2);
        
        Turb_NIR_box_mean = squeeze(nanmean(Turb_NIRsample,[1 2]))';
        lat_box_mean      = squeeze(nanmean(latsample,[1 2]))';
        lon_box_mean      = squeeze(nanmean(lonsample,[1 2]))';
        
        Turb_red_box_std = std((squeeze(std(Turb_redsample,1,'omitnan'))),1,'omitnan');
        lat_box_std      = std((squeeze(std(latsample,1,'omitnan'))),1,'omitnan');
        lon_box_std      = std((squeeze(std(lonsample,1,'omitnan'))),1,'omitnan');
        
        % TABELAO
        
        tabledata    = array2table([round(Turb_NIR_box_mean',2), round(Turb_NIR_box_std',2),...
                                    round(lon_box_mean',3), round(lat_box_mean',3)], ...
                                    'VariableNames',{'Turb_NIR','Turb_NIR_std','lon_sample','lat_sample'});
        date         = datestr(satdate);
        date         = table(datetime(date,'Format','dd.MM.yyyy'),'VariableNames',{'dd.mm.yyyy'});
        sat_scene    = array2table(string(x2), 'VariableNames',{'Scene'});
        big_table_Landsat_RS4 = [big_table_Landsat_RS4; date, sat_scene,tabledata];
    end

disp(i)

end

big_table_Landsat_RS4(isnan(big_table_Landsat_RS4.Turb_nir),:) =[];
save('Turb_L5L7L8L9_histogram.mat','big_table_Landsat_RS4')


