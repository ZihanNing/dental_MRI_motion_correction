
function [tileSize] = estimateTileSize(k, NSweeps)

%ESIMATETILESIZ estimates the tile size of a DISORDER sampling trajectory.
%   [TILESIZE]=ESIMATETILESIZ(K)
%   * K is the Phase Encoding (PE) sampling trajectory of each readout sample as an matrix of size [NumReadoutSamples NumPEDim].
%   * {NSWEEPS} is the number of sweeps of the acquisition. This can help disambiguate the tile size.
%   ** TILESIZE is the tile size as an array of size [1 NumPEDim].
%
%   Yannick Brackenier 2023-03-23

if nargin<2 || isempty(NSweeps);NSweeps=[];end

%%% COMPUTE THE SAMPLING DIFFERENCE
dk = diff(k,1);
nDim = size(k,2);

%%% FIND TILE SIZE
tileSize = zerosL(dk(1,:));

for dim=1:nDim
    %Find unique values
    [dkUnique] = unique(dk(:,dim));
    
    elementsFound = zerosL(dkUnique);
    for ii = 1:length(dkUnique)
        idx = dk(:,dim)==dkUnique(ii);
        elementsFound(ii) = sum(idx);
    end
    assert(sum(elementsFound)==size(dk,1),'Unique values not consistent.');

    %Extract largest dk as the tile size
    [~,idxSort] = sort(elementsFound,'descend');
    tileSize(dim) = abs(dkUnique(idxSort(1)));

    %Check that tile size is not zero for in case sequential sampling is present in (part of) one PE dimension
    if tileSize(dim)==0
        tileSize(dim) = abs(dkUnique(idxSort(2)));
    end
end

%%% CHECK IF IT IS CONSISTENT WITH NUMBER OF SWEEPS PROVIDED
if ~isempty(NSweeps)
    if prod(tileSize)~=NSweeps
        [~,dimWrong] = min(tileSize);
        dimCorrect = dynInd(1:nDim, ~ismember(1:nDim,dimWrong),2);
        tileSize(dimWrong) = NSweeps/tileSize(dimCorrect);
    end
end