
function [co] = defaultColors(n, makeCell)

%%% DEFAULTCOLORS creates a matrix of colors to be used in plots.
%
%   [CO] = DEFAULTCOLORS(N)
%   * {N} is the number of colors to output.
%   * {MAKECELL} is a flag to store all the colors in a cell array.
%   ** CO are the colors as a matrix of size [N 3] or a cell array each element containing a color 
%
%   Yannick Brackenier 2022-11-30

if nargin<1 || isempty(n);n=[];end
if nargin<2 || isempty(makeCell);makeCell=0;end

%%% Manually defined color set
co=[0.8500    0.3250    0.0980;         
    0         0.4470    0.7410;    
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840
    0.4470    0.7410    0
    0.7410    0         0.4470
    0         0         0];

%%% Extract
if ~isempty(n)
    assert(n<=size(co,1),'defaultColors:: Number of requested colors bigger than pre-defined default values. Manually add more here.');
    co = dynInd(co,1:n,1);
end

%%% Make cell
if makeCell
    coNew = [];
    for i=1:size(co,1);coNew{i}=co(i,:);end
    co=coNew;
    coNew=[];
end

