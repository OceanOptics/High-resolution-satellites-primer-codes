%-------------------------------------------------------------------------%
%                                                                         %
%  This code estimates SPM + uncentinty from Water Reflectance at multiple
%   wavelenghts (bands) using the MW algorithm (Tavora et al., 2020)      %
%
%                                                                         %
%                   --------------------------------                      %
%   You need:
%  1. Landsat 8 scene (USGS, level 2 collection 2)                        %
%  2. the MW algorithm package 
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


pathDir = '/pathTo/svalbard_files/usgs_l2c2/'; %change directory to where usgs level2 collection 2 file is
list_scenes = dir([pathDir 'L*']);


match_ups = [];
for ii=1:length(list_scenes)
    
    fileName = char({list_scenes(ii).name});
    imageFile = [pathDir fileName];

    % Run the function
    [dt, lat, lon, RRS, fmask] = getLandsatL2_RW(imageFile, false, [0 99]);
    
    
    [x,l] = find( (year(dt(:,1)) == year(data.date)) & ...
        (month(dt(:,1)) == month(data.date)) & ...
        (day(dt(:,1)) == day(data.date)));
    
    
    field_time = data.time(x); field_time.Format = 'HH:mm:ss';
    
    time_dif = abs((dt.Hour.*60 + dt.Minute) - ...
        (field_time.Hour.*60 + field_time.Minute));
    
    data_select = data(x,:);
    
    for k = 1:size(x,1)
        [row(k), column(k), ~] = findClosestPixel(lon, lat, data_select.lon(k), data_select.lat(k));
        RW_select(k,:) = squeeze(RRS(row(k),column(k),:));
        
        
        temp = data_select.SST; id = isnan(temp); temp_median = median(temp(~id)); temp(id) = repmat(temp_median,size(temp(id),1),1);
        
        [U, rrs, nm] = RW2rrs(RW_select(k,:), wavelengths, L);
        [IOP_matrix,a_sw,dim] = gen_IOP_array(rrs,U,rrs.*NaN,bandwidth,SNR,nm,temp,S,Y,ap443, ap750,bbp700);
        [SPM_wm(k),err_wm(k)] = MW_algorithm(nm, rrs.*NaN, rrs, U, IOP_matrix, a_sw, Q_filter, SNR, dim, temp);
        
    end
    scatter(data_select.SPM,SPM_wm)
    hold on
    
    
    %scene limits
    if ii==1
        
        square_RW             = [reshape(RRS(3329:3620,2563:2855,1),[],1),...
                                 reshape(RRS(3329:3620,2563:2855,2),[],1),...
                                 reshape(RRS(3329:3620,2563:2855,3),[],1)]; %reducing the size of scene to area of interest
       
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
    date         = table(repmat(datetime(dt,'Format','HH:mm:ss'),size(x,1),1),...
        'VariableNames',{'sat_HH.mm.ss_UTC'});
    
    match_ups = [match_ups; date, data_select, tabledata];
    
    clear data_select SPM_wm err_wm
end

%-------------------------------------------------------------------------%




