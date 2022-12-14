function [Rrs_conv, central_band] = conv2sat(SRF_centralband,SRF_curve,Lw, Ed)

%% prepare in-situ radiometric data

[nr_bnds_in,m,radiometric_insitu_Lw] = swapDim(Lw);
[~,~,radiometric_insitu_Ed]          = swapDim(Ed);

nm_insitu              = radiometric_insitu_Lw(:,1);
radiometric_insitu_Lw  = radiometric_insitu_Lw(:,2:end);
radiometric_insitu_Ed  = radiometric_insitu_Ed(:,2:end);

nrvar                  = m-1; clear m

sband               = nm_insitu(1);
eband               = nm_insitu(nr_bnds_in);
nm_insitu_interp    = sband:1:eband;

radiometric_insitu_interp_Lw = interp1(nm_insitu,radiometric_insitu_Lw,nm_insitu_interp);
radiometric_insitu_interp_Ed = interp1(nm_insitu,radiometric_insitu_Ed,nm_insitu_interp); clear radiometric_insitu_Ed radiometric_insitu_Lw

%% prepare satellite related data

% spectral response function of satellite sensor

BNDnm_insitu_interSNR_srf  = load(SRF_centralband);
SRF0                       = load(SRF_curve);

central_band              = round(BNDnm_insitu_interSNR_srf(:,1)); clear BNDnm_insitu_interSNR_srf SRF_center SRF_file
[nn,~,SRF0]               = swapDim(SRF0);
nm_SRF                    = SRF0(:,1);


%% prepare for convolution

SameBands   = ismember(central_band,nm_insitu_interp);
ibnd        = find(SameBands);
nr_bnds_SRF = length(ibnd);

if  nr_bnds_SRF==0
    DISP("no cooresponding bands!!")
    return
end


bnds_SRF   = central_band(SameBands);
ix1        = find(nm_SRF==sband);
ix2        = find(nm_SRF==eband);
nm_SRF     = nm_SRF(ix1:ix2);

SRF0       = SRF0(ix1:ix2,ibnd+1);



for i=1:nr_bnds_SRF
    
    index=find(SRF0(:,i)>0);
    
    if ~isempty(index)
        
        nm_conv=nm_SRF(index);
        
        SRF   =  SRF0(index,i);
        radiometric_insitu_data_Lw  =  interp1(nm_insitu_interp,radiometric_insitu_interp_Lw,nm_conv);
        radiometric_insitu_data_Ed  =  interp1(nm_insitu_interp,radiometric_insitu_interp_Ed,nm_conv);

        radiometric_insitu_conv_Lw(:,i)  =  (radiometric_insitu_data_Lw'*SRF)./nansum(SRF);
        radiometric_insitu_conv_Ed(:,i)  =  (radiometric_insitu_data_Ed'*SRF)./nansum(SRF);

        
    end
end

Rrs_conv = Lw./Ed; 


end

function [n,m,radiometric_insitu]=swapDim(radiometric_insitu)

[n m] =size(radiometric_insitu);

if n<m
    radiometric_insitu=radiometric_insitu';
    swap=m;
    m=n;
    n=swap;
    
end

end

