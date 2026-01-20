
function [binSample] = makeBins(N, NBins)

%MAKEBIN creates an index array linking samples to equidistant bins.
%   [BINSAMPLE]=MAKEBIN(N,{NBINS})
%   * N the number of samples of the array. If an array, its number of samples will be derived from it.
%   * {NBINS} the number of bins to use to generate the index list.
%   ** BINSAMPLE is the index list.
%
%   Yannick Brackenier 2023-07-05

if nargin<2 || isempty(NBins);NBins = length(N);end

if numel(N)>1;N = length(N);end

binSample=ceil(NBins*(((1:N)-0.5)/N));
