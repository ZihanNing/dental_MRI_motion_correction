
function [x, cMap] = mapToColormap(x, cmapInfo, rangeAssign, numColors)

%MAPTOCOLORMAP maps numerical values in an array to a specified colormap.
%   [X,CMAP]=MAPTOCOLORMAP(X, CMAPINFO, {RANGEASSIGN},{NUMCOLORS})
%   * X is the array to assign colors to.
%   * CMAPINFO is the colormap information. Can be scalar/string/matrix/colormap (see below). Defaults to jet.
%   * {RANGEASSIGN} is the specific range of values in x to map onto the  specifiec colormap.
%   * {NUMCOLORS} is the number of colors to use in the colormap in case CMAPINFO is a scalar.
%   ** X is the (flattened??) array number of elements that belong to any interval.
%
%   Yannick Brackenier 2022-08-30

if nargin<2 || isempty(cmapInfo); cmapInfo = 1;end
if nargin<3 || isempty(rangeAssign); rangeAssign = [min(x(~isinf(x)))  max(x(~isinf(x)))];end
if nargin<4 || isempty(numColors); numColors = 500;end

N = size(x);
N = N(1:numDims(x));

%%% Create colormap - TODO: move to a separate function
if isscalar(cmapInfo)
    if cmapInfo==0;cMap = gray(numColors);
    elseif cmapInfo== 1; cMap = jet(numColors);
    elseif cmapInfo== 2; cMap = turbo(numColors);%Introduced in 2020
    elseif cmapInfo== 3; cMap = hot(numColors);
    elseif cmapInfo== 4; cMap = cmrmap(numColors);
    else; error('mapToColormap:: cmapInfo scalar not implemented.'); 
    end
end
if ismatrix(cmapInfo) && size(cmapInfo,1)>1; cMap = cmapInfo; assert(size(cMap,2)==3,'mapToColormap:: cmapInfo matrix containing colormap not valid.');end
if ischar(cmapInfo); cMap = colormap(cmapInfo);end
if iscell(cmapInfo) ; cMap = createColormap(cmapInfo,numColors);end

%%% Rescale array from a range to the size of the colormap
x = rescaleND(x,[1 size(cMap,1) ], rangeAssign);
x = max(x,1);
x = min(x,size(cMap,1));

%%% Assign colors to samples
x = x(:);%Flatten
x = cMap(round(x),:);%Assign color
x = resSub(x,1,N);%Reshape back so that there are now N+1 dimensions with the last one the RGB one
