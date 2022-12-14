function [dt, lat, lon, RW_stack, fmask] = getLandsatL2_RW(pathfolder, retrieve_land, prc_lim)
% code adapted by Juliana Tavora in November 2022
% original code for SST by Guillaume Bourdin
% created: April 12, 2021
%%
% Extract Water reflectance from Landsat 5, 7 and 8 collection 2 level 2 images and apply masks
%
% INPUT:
%   - pathfolder: Landsat-8 folder path <1xM char>
% Optional input:
%   - retrieve_land: <boolean> retrieve land; default = false
%   - prc_lim: <1x2 double> percentils to remove (%); default = [2.5 99]
%
% OUTPUT:
%   - dt: <1x1 datetime> satellite over pass date and time
%   - lat: <NxM double> matrix of latitudes
%   - lon: <NxM double> matrix of longitudes
%   - RW_stack: <NxM double> matrix of RW
%   - fmask: <NxM double> matrix of masks generated from fmask
%
% examples:
%    - [dt, lat, lon, temperature, fmask] = getLandsatL2_RW(pathfolder)
%    - [dt, lat, lon, temperature, fmask] = getLandsatL2_RW(pathfolder, true)
%    - [dt, lat, lon, temperature, fmask] = getLandsatL2_RW(pathfolder, false, [5 97.5])
%%
if nargin < 1
  error('Not enough input argument')
elseif nargin < 2
  fprintf('Warning: Land, small lakes and rivers retrieval not defined, default = false\n')
  retrieve_land = false;
  fprintf('Warning: Percentils limit not defined, default = [2.5 99]\n')
  prc_lim = [2.5 99];
elseif nargin < 3
  fprintf('Warning: Percentils limit not defined, default = [2.5 99]\n')
  prc_lim = [2.5 99];
elseif nargin > 3
  error('Too many input arguments,')
end

% disable warnings
warning('off')

% list files in directory
[~, filename] = fileparts(pathfolder);
foo = {dir(fullfile(pathfolder, 'L*')).name}'; % (1:end-4)
foo = cellfun(@(x) [pathfolder filesep x], foo, 'un', 0); % (1:end-4)

% read MTL file
MTLfile = fileread([foo{contains(foo, '_MTL.txt')}]);

% get variables from MTL files
spacecraft = parseMTL(MTLfile, 'SPACECRAFT_ID =');

% set thermal band by spacecraft (B6 for L4, L5, L7 and B10 for L8)
switch spacecraft
  case {'LANDSAT_4','LANDSAT_5','LANDSAT_7'}
    STBAND = '6';
  case {'LANDSAT_8','LANDSAT_9'}
    STBAND = {'4','5','6'};
  otherwise
    error('Spacecraft %s not known', spacecraft)
end

% get calibration slope and intercept
for i=1:size(STBAND,2)
    slope_char = parseMTL(MTLfile, ['REFLECTANCE_MULT_BAND_' STBAND{i} ' =']); slope(i) = str2num(slope_char(1:8));
    intercept_char = parseMTL(MTLfile, ['REFLECTANCE_ADD_BAND_' STBAND{i} ' =']); intercept(i) = str2num(intercept_char(1:4));
end
clear intercept_char slope_char i
SUN_EL = parseMTL(MTLfile, 'SUN_ELEVATION ='); 

% extract date time and info for lat/lon
date = parseMTL(MTLfile, 'DATE_ACQUIRED =');
time = parseMTL(MTLfile, 'SCENE_CENTER_TIME =');
dt = datetime([date ' ' time], 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSSS''Z');



if retrieve_land
  fprintf('Applying ice and cloud masks to %s ... ', filename)
else
  fprintf('Applying land, ice and cloud masks to %s ... ', filename)
end
% reconstruct fmask output Quality Assessment Band
fmask = qa_to_mask(foo{contains(foo, '_QA_PIXEL.TIF')});
fprintf('Done\n')


% getband data
fprintf('Loading band from %s ... ', filename)
for i = 1:size(STBAND,2)
    band = [foo{contains(foo, ['_SR_B' STBAND{i} '.TIF'])}];
    t = Tiff(band);
    RW = double(read(t));
    RW_corr = (RW * slope(i) + intercept(i))./sin(SUN_EL);
    
    % apply fmask
    if retrieve_land
        RW_corr(fmask > 1) = NaN;
    else
        RW_corr(fmask ~= 1) = NaN;
    end
    
    RW_stack(:,:,i) = RW_corr;
end
fprintf('Done\n')


% Remove aberrant data
RW_stack(RW_stack < 0) = NaN;

% get image lat/lon
fprintf('Recovering latitude and longitude from %s ... ', filename)
info = geotiffinfo(band);
img_height = info.Height; % Integer indicating the height of the image in pixels
img_width = info.Width; % Integer indicating the width of the image in pixels
[cols,rows] = meshgrid(1:img_width, 1:img_height);
[x,y] = pix2map(info.RefMatrix, rows, cols);
[lat,lon] = projinv(info, x, y);
fprintf('Done\n')




% Remove data out of predefined percentils (default = [2.5 99])
prct = prctile(RW_stack(:), prc_lim);
fprintf('%.2f%% pixel deleted from %s: out of [%.1f %.1f] percentils\n', ...
  sum(sum(RW_stack < prct(1) | RW_stack > prct(2))) / size(RW_stack(:), 1) * 100, ...
  filename, prc_lim(1), prc_lim(2))
RW_stack(RW_stack < prct(1) | RW_stack > prct(2)) = NaN;
end


function parsed_parameters = parseMTL(MTLfile, toparse)
% Author: Guillaume Bourdin
% Date: Jan 20th, 2021
%
% Parse MTL file of landsat
if ~iscell(toparse)
  parsed_parameters = strrep(cell2mat(regexp(MTLfile, ['(?<=' toparse ' ).(.*)'],...
    'match','dotexceptnewline')), '"', '');
  try
    foo = str2double(parsed_parameters);
    if ~isnan(foo)
      parsed_parameters = foo;
    end
  catch
  end
else
  parsed_parameters = cell(size(toparse));
  for i = 1:max(size(toparse))
    parsed_parameters{i} = strrep(regexp(MTLfile, ['(?<=' toparse{i} ' ).(.*)'],...
      'match','dotexceptnewline'), '"', '');
    try
      foo = str2double(parsed_parameters{i});
      if ~isnan(foo)
        parsed_parameters{i} = foo;
      end
    catch
    end
  end
end
end

%% function qa_to_mask to generate f_mask output from QA band
% Author: Gabriel Hesketh
% Date: April 2021

% A function to read a level 2, collection 2 Landsat-8 file's QA band
% and generate a working mask for the file.
% The only band necessary for mask generation is qa_pixel.TIF; however,
% creating images from this file also requires the accompanying MTL file.
% qa_path          The qa_pixel band's full path 
% return f_mask   The generated mask for the L2C2 file
function f_mask = qa_to_mask(qa_path)
% here is an example qa_path, comment out the function line and corresponding end 
% if you want to do make your own list or processes scenes individually

info = Tiff(qa_path, 'r'); %read the file
qa = read(info);
f_mask = zeros(size(qa, 1), size(qa, 2)); %make a blank array to fill

% now get each bit of the band to get the important values from the band
bit1 = bitget(qa, 1); %fill
bit2 = bitget(qa, 2); %dilated cloud
bit3 = bitget(qa, 3); %cirrus
bit4 = bitget(qa, 4); %cloud 
bit5 = bitget(qa, 5); %cloud shadow
bit6 = bitget(qa, 6); %snow
% bit7 = bitget(qa, 7); %clear
bit8 = bitget(qa, 8); %water

% fill the zeros array based on the bit bands

% from FMask 4.3's autoFmask file, this shows the values given
% autoFmask, 
%     fmask      0: clear land
%                1: clear water
%                2: cloud shadow
%                3: snow
%                4: cloud
%                5: filled (outside)*

% (in autoFmask, outside values are set to 255)

f_mask(bit8 == 1) = 1; %set water
f_mask(bit2 == 1) = 4; %set clouds
f_mask(bit3 == 1) = 4; 
f_mask(bit4 == 1) = 4; 
% f_mask(bit7 == 0) = 4; 
f_mask(bit5 == 1) = 2; %overrwite some of cloud with cloud shadow
f_mask(bit6 == 1) = 3; %set snow
f_mask(bit1 == 1) = 5; %set outside
end
