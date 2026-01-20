

function [filterSize, isSlab, dimSlab] = filterForSlab(FOVmm, factor)

%FILTERFORSLAB detect the slab dimension from a FOV and returns a filter used for slab apodisation during motion correction.
%   [FILTERSIZE,ISSLAB,IDX]=FILTERFORSLAB(FOV,{FACTOR})
%   * FOV is the FOV of the array in a.u..
%   * {FACTOR} is the factor of the maximum FOV dimension to detect the slab dimension (if it's lower).
%   ** FILTERSIZE is the normalised cutoff of the filter (normalised to 1).
%   ** ISSLAB is a flag to indicate a slab was detected.
%   ** IDX is the dimension of the slab.
%
%   Yannick Brackenier 2023-03-17

if nargin<2 || isempty(factor);factor=.4;end

%%% DETECT SLAB
isLowFOV = FOVmm<factor*multDimMax(FOVmm);
isSlab = any(isLowFOV);

%%% SET FILTER DIMENSION
if isSlab 
    if multDimSum(isLowFOV)>1
        warning('Slab in 2 direction is not expected and hence filter is disabled.');
        dimSlab = [];
        isSlab = 0;
    else
        dimSlab = find(isLowFOV==1);
    end   
else
    dimSlab = [];
end

filterSize = [20 20 20];
filterSize(dimSlab)=1;
