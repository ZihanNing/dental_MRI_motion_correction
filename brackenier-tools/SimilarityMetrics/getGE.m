
function grEn=getGE(x,ROI)

%GRADIENTENTROPY  Computes the gradient entropy of a multidimensional
%array
%   GREN=GRADIENTENTROPY(X)
%   * X is the multdimensional array
%   ** GREN is the gradient entropy
%

if nargin<2 || isempty(ROI);ROI=[];end

ND=numDims(x);
for n=1:ND
    xd{n}=abs(x-dynInd(x,[2:size(x,n) 1],n));
    xd{n}=xd{n}+1e-12;
end
xd=cat(ND+1,xd{:});

%%% WEIGHT BY ROI
if ~isempty(ROI); xd=xd.*single(ROI);end    

xd=xd/(sum(xd(:)));  

grEn = sum(xd(:).*log2(xd(:)));%Maximum possible entropy, all outcomes equally probable
if  isempty(ROI)
    grEn = grEn/log2(numel(xd));
else
    grEn = grEn/log2(multDimSum(ROI));
end
grEn = gather(grEn);