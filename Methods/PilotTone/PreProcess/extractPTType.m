

function p = extractPTType(p, useRealImag, NChaRec)

if nargin<3 || isempty(NChaRec);NChaRec=[];end
isSlice = size(p,4)>1;
if isSlice; dimCha=4;dimSamples=1:2; else; dimCha=1;dimSamples=2;end

%%% Class of signal
if useRealImag == 1
    p = cat(dimCha, real(p), imag(p));
end

%%% Number of channels
if ~isempty(NChaRec)
    p = dynInd(p,1:NChaRec, dimCha);
end


