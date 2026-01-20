
function [xM, coverage, shift] = recenterImage(x, dim, extractCoverage, M)

%RECENTERIMAGE   re-centers an image to the logical center of the array.
%   [X,COVERAGE,SHIFT]=RECENTERIMAGE(X, DIM, {EXTRACTCOVERAGE})
%   * X the array to recenter.
%   * DIM the dimensions of the array to recenter.
%   * {EXTRACTCOVERAGE} extracts the part of the array that only contains non-non-zero elements.
%   * {M} the mask to use for centroid computation instead of a default-computed one.
%   ** XM is the recentered array.
%   ** COVERAGE the range in that dimensions that has non-zero values.
%   ** SHIFT the shift that was applied.
%
%   Yannick Brackenier 2023-04-17

if nargin<3 || isempty(extractCoverage);extractCoverage=0;end
if nargin<4 || isempty(M);M=[];end

if numDims(x)>3
    xM = multDimMea(abs(x),4:numDims(x));
    warning('recenterImage:: Dimensions >4 averaged to compute center on single 3D volume. TO DO: extend to each volume on ND array.')
else
    xM = abs(x);
end

%%% General parameters
deb = 0;
if deb;xOrig=xM;end
N = size(xM);
nD = numDims(xM);

%%% Create mask
if isempty(M)
    parS=[];
    parS.Otsu = 0:.25:1;
    parS.nDilate = 0;
    M = refineMask(xM,parS);
end

%%% Average dimensions not to recenter
dimAverage = 1:nD;
dimAverage(ismember(dimAverage,dim)) = [];
if ~isempty(dimAverage);M = multDimSum(M,dimAverage);end

%%% Find centroid in each dimension
coverage = zeros([nD 2]);
center = zeros([nD 1]);

for i=1:length(dim)
    %Average
    dimAverage = 1:nD;
    dimAverage(ismember(dimAverage,dim(i))) = [];
    m = multDimMea(M,dimAverage);
    %Find center
    [~,center(i)] = max(m(:));
    %Refine center by Gaussian fit on peak
    mTemp = dynInd(m(:), (center(i)-round(N(dim(i))/4)):(center(i)+round(N(dim(i))/4)),1);
    [mTempFit, ~, mu] = gaussianFit(mTemp(:).', 1:length(mTemp),2);
    if deb>1
        figure; hold on; 
        plot(mTemp);
        plot(mTempFit);
    end
    if ~isnan(mu); center(i) = center(i) + round(mu-round(N(dim(i))/4));end
    %Find coverage
    tt = find(m(:));
    coverage(i,:) = [tt(1) tt(end)];
end

%%% Apply recentering
shift = zeros([1 nD]);
shift(dim) = center(dim).' - centerIdx(N(dim));
x = circshift(x,-shift);
if deb
    xM = circshift(xM,-shift);
    plotND([],cat(4,xOrig,xM),[],[],0,[],[],{'Original';'Centered'},onesL(xM),{2});
end

%%% Adjust coverage and make correction
coverage = coverage - repmat(shift(:),[1 2]);
for i=1:length(dim)
   coverage(i,1) = max(coverage(i,1),1);
   coverage(i,2) = min(coverage(i,2),size(xM,dim(i)));
end

%%% Extract
if extractCoverage
    for i=1:length(dim); xM = dynInd(xM, coverage(i,1):coverage(i,2) ,dim(i)); end
end
