
function [x,idx] = fillCent (x, dim, val, fillLower, fillUpper)

%FILLCENT   fills the center (in the DFT sense) of an array.
%   [X, IDX] = FILLCENT(X,{DIM},{VAL},{INCLUDELOWER},{INCLUDEUPPER})
%   * X is the array
%   * {DIM} is the set of dimensions along which to fill the array
%   * {VAL} is the value to fill with
%   * {INCLUDELOWER} whether also to fill the lower 1/4 of the array
%   * {INCLUDEUPPER} whether also to fill the upper 1/4 of the array
%   ** X is the filled array
%   ** IDX is a cell with each element containing the indices of the array along that dimension that has been filled
%
%   Yannick Brackenier 2022-07-17

if nargin<2 || isempty(dim);dim=1:numDims(x);end
if nargin<3 || isempty(val);val=0;end
if nargin<4 || isempty(fillLower);fillLower=0;end
if nargin<5 || isempty(fillUpper);fillUpper=0;end

N=size(x);
zeroF=ceil((N+1)/2);
idx=[];

%%% Compute support
for i=1:length(dim)
    %Default center indices
    orig=zeroF(dim(i))-ceil(((N(dim(i))/2)-1)/2);
    fina=zeroF(dim(i))+floor(((N(dim(i))/2)-1)/2);
    %Include excpetions
    if fillLower;orig(1)=1;end
    if fillUpper;fina(1)=N(dim(i));end
    %Create indices
    idx{i} = orig(1):fina(1);
end

%%% Fill
if ~isempty(dim)
    x=dynInd(x,idx,dim,val);
end
