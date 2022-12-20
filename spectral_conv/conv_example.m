clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example code to spectral conv. in-situ hyperspectral radiometric data
%       * same procedure can be carried out to spectrally conv other
%         hyperspectral data (e.g., algorithm coeffcients)
%
% This piece of code requires: 1. central bands from satellite sensor
%                              2. SRF from satellite sensor
%                              3. in-situ radiometric Lw and Ed data and wavelenghts
%
% Outputs: Rrs spectrally conv. 
%
%
% In-situ data used here is available from the SeaSWIR dataset:  https://doi.pangaea.de/10.1594/PANGAEA.886287
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%seaswir dataset:
Lw = table2array(import_data("/path_to_/datasets/SeaSWIR_TRIOS_Lw.tab", [239, Inf]));
Ed = table2array(import_data("/path_to_/datasets/SeaSWIR_TRIOS_Ed.tab", [239, Inf]));
nm = 350:2.5:900;

Lw = [nm; Lw]'; 
Ed = [nm; Ed]'; clear nm


%% prepare in-situ radiometric data


SRF_centralband  = [442.7 492.7 559.8 664.6 704.1 740.5 782.8 832.8 864.7 945.1 1373.5 1613.7 2202.4];
SRF_curve        = 'MSI_A_SRF.m';



[Rrs_conv, central_band] = conv2sat(SRF_centralband,SRF_curve,Lw, Ed);

plot(nm_insitu,Lw(:,2)./Ed(:,2))
hold on
plot(central_band(SameBands),Rrs_conv(1,:),'o')


