function [ BWAI ] = BWAI_Landsat( Grn, Red, NIR, Lambda_Grn, Lambda_Red, Lambda_NIR, Blue, SWIR1, Lambda_Blue, Lambda_SWIR1, MaskValid )
%%% This code for BWAI calculation was developed by Dr. Yongquan Zhao, E-mail: yqzhao@link.cuhk.edu.hk
%%% Reference:
%%% Zhao, Y., Liu, D., & Wei, X. (2020). Monitoring cyanobacterial harmful algal blooms at high spatiotemporal resolution by fusing Landsat and MODIS imagery. Environmental Advances, 2, 100008. doi:10.1016/j.envadv.2020.100008
%%% Last modification date: December 10, 2020.

% Inputs:
% Blue:  the Landsat blue band.
% Grn:   the Landsat green band (the reflectance peak of chlorohpyll (subsurface algae)).
% Red:   the Landsat red band.
% NIR:   the Landsat near-infrared band (the reflectance peak of floating algae).
% SWIR1: the Landsat shortwave infrared band 1.
% Lambda_Blue:  the central wavelength of the Landsat blue band.
% Lambda_Grn:   the central wavelength of the Landsat green band.
% Lambda_Red:   the central wavelength of the Landsat red band.
% Lambda_NIR:   the central wavelength of the Landsat NIR band.
% Lambda_SWIR1: the central wavelength of the Landsat SWIR1 band.
% MaskValid: the mask for valid pixels (clouds and dead/noise pixels are masked).
% PS: All the input water surface reflectance bands are atmospheric corrected and multiplied with a scale factor of 10000.

[H, W] = size(Grn);

% Detecting green or NIR reflectance peaks.
RPH = zeros(H, W);
Lambda = [Lambda_Grn, Lambda_NIR];
Th = 600;
for i=1:H
    for j=1:W
        if (MaskValid(i,j) == 0)
            continue;
        end
        
        if(Blue(i,j)>Th) % Exclude false positives initially.
            RPH(i,j) = NaN;
        else
            Rho = [Grn(i,j), NIR(i,j)];
            [RhoMax, MaxInd] = max(Rho);
            Lambda_Max = Lambda(MaxInd);
            RPH(i,j) = RhoMax - (Blue(i,j) + (SWIR1(i,j)-Blue(i,j)) * ((Lambda_Max-Lambda_Blue)/(Lambda_SWIR1-Lambda_Blue)));
        end
    end
end

% The RPH depressing factor (for high reflectance at green bands due to water suspended sediments).
Fs = Red - (Grn + (NIR-Grn) * ((Lambda_Red-Lambda_Grn)/(Lambda_NIR-Lambda_Grn)));

% The mask for pixels with high suspended sediment loadings, 30 is a empirical threshold (the theoretical value is 0).
T = 30;
SediMask = (Fs>T);

% The RPH modulating factor (for waters with CyanoHABs and low suspended sediment loadings).
Fc = CalNDXI( Grn, Blue, MaskValid );

% Depressing and modulating RPH.
BWAI = (RPH .* exp(Fc) .* (1-SediMask) + RPH ./ exp(Fs/10000) .* SediMask) .* MaskValid;



function NDXI = CalNDXI( B1, B2, MaskValid )
% The range of NDXI is [-1, 1].

% Offset negative values in Landsat bands.
MinVal1 = min(min(B1.*MaskValid));
MinVal2 = min(min(B2.*MaskValid));
MinVal = min([MinVal1 MinVal2]);
if MinVal<0
    B1 = B1 - MinVal;
    B2 = B2 - MinVal;
end

NDXI = (B1 - B2) ./ (B1 + B2); % NDXI.
return;

