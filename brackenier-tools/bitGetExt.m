
function [bit] = bitGetExt(x, id, base)

%Extension of bitget for arbitrary base
%Yannick Brackenier

assert(mod(base,1)==0,'bitGetExt:: Base must be integer.');

nD=numDims(x);
if nD==0;nD=1;end

maxId = max(id);
bit = [];
for i=1:maxId
    bitTemp = mod(x, base);
    x = (x-bitTemp)/base;
    bit = cat(nD+1,bit, bitTemp);
end

%Extract queried ones
bit = dynInd(bit,id,numDims(bit));