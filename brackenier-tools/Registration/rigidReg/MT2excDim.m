
function [dimFH, fac] = MT2excDim(MT, fac)

if nargin<2 || isempty(fac);fac=.3;end

[permTemp,flTemp] = MT2perm(MT);
dimFH = permTemp(3);
if flTemp(dimFH)==1
    dimFH=-dimFH;
end
