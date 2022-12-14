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
% Outputs: Rrs spectrally conv. for Sentinel 2 A bands, to apply for Sentinel 2B change SRF and central bands
%
%
% In-situ data used here is available from the SeaSWIR dataset:  https://doi.pangaea.de/10.1594/PANGAEA.886287
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%seaswir dataset:
Lw = table2array(import_data("/Users/julianatavora/Documents/convolvi_github/datasets/SeaSWIR_TRIOS_Lw.tab", [239, Inf]));
Ed = table2array(import_data("/Users/julianatavora/Documents/convolvi_github/datasets/SeaSWIR_TRIOS_Ed.tab", [239, Inf]));
nm = 350:2.5:900;

Lw = [nm; Lw]'; 
Ed = [nm; Ed]'; clear nm


%% prepare in-situ radiometric data


SRF_centralband  = 'MSI_centralBands.m';
SRF_curve        = 'MSI_A_SRF.m';



[Rrs_conv, central_band] = conv2sat(SRF_centralband,SRF_curve,Lw, Ed)
