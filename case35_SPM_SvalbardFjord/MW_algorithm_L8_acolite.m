%-------------------------------------------------------------------------%
%                                                                         %
%  This code estimates SPM + uncentinty from Water Reflectance at multiple
%   wavelenghts (bands) using the MW algorithm (Tavora et al., 2020)      %
%
%                                                                         %
%                   --------------------------------                      %
%   You need:
%  1. Landsat 8 scene (acolite processed, level 2 collection 2)           %
%  2. the MW algorithm package: https://github.com/OceanOptics/MW_algorithm_satellite
%  3. in-situ data 'in-situ_data_svaalbard.mat' available in the root folder

% --                                                                      %
% Juliana Tavora, 2019                                                    %
%                                                                         %
%-------------------------------------------------------------------------


close all
clear all
clc

%-------------------------------------------------------------------------%

%set insitu data
load('in-situ_data_svaalbard.mat')
data = [Kronebreen; Tunabreen];


%-------------------------------------------------------------------------%

addpath '/pathTo/MW_algorithm/' %change directory


% inputs necessary for MW algorithm 

% gordon et al (1988)
L = [0.0949,0.0794];
Q_filter = 0.5;

%slopes/coefs used in the inversions
S         = 0.006:0.001:0.014;
Y         = 0:0.1:1;
ap443     = 0.01:0.01:0.06;
ap750     = 0.013:0.001:0.015; 
bbp700    = 0.002:0.001:0.021;

SNR = [227; 201; 267];
bandwidth = [37; 28; 85];
wavelengths = [655 865 1609];

%-------------------------------------------------------------------------%

pathDir = 'pathto/svalbard/acolite_l2c2/';
list_scenes = dir([pathDir '*.nc']);


match_ups = [];
for ii=1:length(list_scenes)
    
    fileName = char({list_scenes(ii).name});
    imageFile = [pathDir fileName];
    
    
    %latitude/long from the image
    lat = ncread(imageFile, 'lat'); lat = lat';
    lon = ncread(imageFile, 'lon'); lon = lon';
    
    mask =  ncread(imageFile, 'rhos_2201');
    
    satdate = datetime(str2double(fileName(8:11)),str2double(fileName(13:14)),str2double(fileName(16:17)),'Format','dd.MM.yyyy');
    
    sat_time = ncreadatt(imageFile,'/','isodate');
    sat_time = datetime(sat_time(12:19), 'Format', 'HH:mm:ss');
    
    [m,n] = size(lat);
    
    %%% Reading in Rrs_vvv into one 3D array, just easier to work with
    RRS = zeros(m,n,length(wavelengths));
    for i = 1:length(wavelengths)
        f_index = ['rhos_' num2str(wavelengths(i))];
        rrs_array = ncread(imageFile, f_index); rrs_array(rrs_array<0)=NaN;
        RRS(:,:,i) = rrs_array';
    end
    
    
    [x,l] = find( (year(satdate(:,1)) == year(data.date)) & ...
                 (month(satdate(:,1)) == month(data.date)) & ...
                   (day(satdate(:,1)) == day(data.date)));
    
    
    field_time = data.time(x); field_time.Format = 'HH:mm:ss';
    
    time_dif = abs((sat_time.Hour.*60 + sat_time.Minute) - ...
        (field_time.Hour.*60 + field_time.Minute));
    
    data_select = data(x,:);
    
    for k = 1:size(x,1)
        [row(k), column(k), ~] = findClosestPixel(lon, lat, data_select.lon(k), data_select.lat(k));
        RW_select(k,:)         = squeeze(RRS(row(k),column(k),:));
        
        
        temp = data_select.SST; id = isnan(temp); temp_median = median(temp(~id)); temp(id) = repmat(temp_median,size(temp(id),1),1);
        
        [U, rrs, nm] = RW2rrs(RW_select(k,:), wavelengths, L);
        [IOP_matrix,a_sw,dim] = gen_IOP_array(rrs,U,rrs.*NaN,bandwidth,SNR,nm,temp,S,Y,ap443, ap750,bbp700);
        [SPM_wm(k),err_wm(k)] = MW_algorithm(nm, rrs.*NaN, rrs, U, IOP_matrix, a_sw, Q_filter, SNR, dim, temp);
        
    end
    scatter(data_select.SPM,SPM_wm)
    hold on
    
    if ii==3
        
        square_RW             = [reshape(RRS(3329:3620,2563:2855,1),[],1),...
                                reshape(RRS(3329:3620,2563:2855,2),[],1),...
                                reshape(RRS(3329:3620,2563:2855,3),[],1)];
       
         clear RRS U rrs RW_sel* IOP_matrix*
                      
        lat_lala              = lat(3329:3620,2563:2855);
        lon_lala              = lon(3329:3620,2563:2855);
        
        SST                   = nanmedian(temp);
        
        [U, rrs, nm]          = RW2rrs(square_RW, wavelengths, L);
        [IOP_matrix,a_sw,dim] = gen_IOP_array(rrs,U,NaN(1,3),bandwidth,SNR,nm,SST,S,Y,ap443, ap750,bbp700);
        
        SPM_wm_day1           = NaN(size(square_RW,1),1);
        err_day1              = NaN(size(square_RW,1),1);
        
        tic
        for h = 1:size(square_RW,1)
            
            [SPM_wm_day1(h),err_day1(h)] = MW_algorithm(nm, NaN(1,3), rrs(h,:),...
                U(h,:), IOP_matrix, a_sw, Q_filter, SNR, dim, SST);
        end
        toc
         SPM_alg = reshape(SPM_wm_day1,292,293);     SPM_alg_err = reshape(err_day1,292,293);   
    end
            
    
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
    %                        create table with data                           %
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
    
    tabledata    = array2table([round(SPM_wm',2),round(err_wm',2)],...
        'VariableNames',{'SPM_mw','SPM_std'});
    date         = table(repmat(datetime(sat_time,'Format','HH:mm:ss'),size(x,1),1),...
        'VariableNames',{'sat_HH.mm.ss_UTC'});
    
    match_ups = [match_ups; date, data_select, tabledata];
    
    clear data_select SPM_wm err_wm
end

