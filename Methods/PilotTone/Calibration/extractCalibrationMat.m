
function [AfOut,AbOut] = extractCalibrationMat(AfFull, AbFull, Af, Ab, NChaPTRec, calibrationOffset, di)

if nargin<7 || isempty(di);di=1;end

if di==1%Extraction
    AfOut = dynInd(AfFull,1:NChaPTRec, 1);
    %if calibrationOffset; AfOut = cat(2,AfOut, dynInd(AfFull,size(AfFull,2),2));end

    AbOut = dynInd(AbFull,1:NChaPTRec, 2);
    if calibrationOffset; AbOut = cat(2,AbOut, dynInd(AbFull,size(AbFull,2),2));end

else %Filling
    AfOut = dynInd(AfFull,1:NChaPTRec, 1, Af);
    %if calibrationOffset; AfOut = dynInd(AfOut, size(AfFull,1),1, dynInd(Af, size(Af,1),1)  );end
    
	AbOut = dynInd(AbFull,1:NChaPTRec, 2, dynInd(Ab,1:NChaPTRec,2) );
    if calibrationOffset; AbOut = dynInd(AbOut, size(AbFull,2),2, dynInd(Ab, size(Ab,2),2)  );end
end