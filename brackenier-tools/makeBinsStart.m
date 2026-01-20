
function [binSample] = makeBinsStart(N, idStart)

%MAKEBIN creates an index array linking samples to equidistant bins.
%   [BINSAMPLE]=MAKEBINSTART(N,IDSTART)
%   * N the number of samples of the array. If an array, its number of samples will be derived from it.
%   * IDSTART the .
%   ** BINSAMPLE is the index list.
%
%   Yannick Brackenier 2023-11-28

if nargin<2 || isempty(idStart);error('makeBinsStart:: Start indices required');end
idStart = idStart(:).';
if idStart(1)>1;idStart = cat(2, 1, idStart);end

if numel(N)>1;N = length(N);end
binSample=zeros([1 N],'single');

NBins = length(idStart);

for i=1:NBins
    if i==NBins
        binSample(idStart(i):end) = i;
    else
        binSample(idStart(i):(idStart(i+1)-1))=i;
    end
end
