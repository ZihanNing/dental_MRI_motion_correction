

function [limPlot] = getTLims(T, MS, addOffset, scalingFac)

%GETTILIMS  computes the upper and lower limits of a motion trace array (e.g. for plotting purposes).
%   [LIMPLOT]=GETTILIMS(T,{MS},{ADDERROR},{SCALINGFAC})
%   * T is the motion trace where the last dimension are the 6 rigid motion parameters. States are stored in a dimension before the last one.
%   * {MS} is the resolution of the corresponding image.
%   * {ADDOFFSET} a flag to set a small error in case limPlot is all zeros.
%   * {SCALINGFAC} a the scaling to use for range calculation from the maximum values.
%   ** LIMPLOT is an array with size [2 2] with the limits on each row for translation and rotation respectively.
%
%   Yannick Brackenier 2022-01-30

if nargin<2 || isempty(MS);MS=[1 1 1];end
if nargin<3 || isempty(addOffset);addOffset=1;end
if nargin<4 || isempty(scalingFac);scalingFac=1.1;end

%%% Extract largest values for rotation and translation
nD = numDims(T);
maxTran = max(MS(:))*multDimMax(abs(dynInd(T,1:3,nD)));
maxRot = convertRotation( multDimMax(abs(dynInd(T,4:6,nD))),'rad','deg');

%%% Construct the limits
limPlot = [-maxTran maxTran; ...
           -maxRot  maxRot];

%%% Scale and add offset
limPlot = scalingFac*limPlot;
if addOffset;limPlot=limPlot+repmat([-1e-3 1e-3],[2,1]);end

%%% Gather for visMotion.m usage
limPlot = gather(limPlot);
