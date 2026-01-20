
function NGS = getNGS(x, ROI, pow)

%GETNGS computes the Normalised Gradient Squared (MGS) of an image. The NGS is based on McGee KP, Manduca A, Felmlee JP, Riederer SJ, Ehman RL. Image metric‐based correction (autocorrection) of motion effects: analysis of image metrics.
%   [NGS] = GETNGS(X,{ROI}) 
%   * X the array for which to compute the NGS.
%   * {ROI} the region of interest (ROI) for which to compute the similarity metric.
%   * {POW} the power of the normalised gradient.
%   ** NGS the calculated NGS
%
%   Yannick Brackenier 2022-07-22

if nargin < 2 ; ROI = []; end
if nargin < 3 ; pow = 2; end
nD=numDims(x);

%%% COMPUTE FINITE DIFFERENCES
G=[];
for i=1:nD; G = cat(nD+1,G,FiniteDiff(x,i));end

%%% WEIGHT BY ROI
if ~isempty(ROI); G=G.*single(ROI);end    

%%% NORMALISE
G = abs(G)./multDimSum(abs(G));

%%% FINAL METRIC
NGS = multDimSum(G.^pow);
NGS = gather(NGS);
