
function [Af, Ab] = createCalibrationMat(NChaPTRec, calibrationOffset, useRealImage)

Af = zeros([NChaPTRec, 6+calibrationOffset        ],'single');%6 motion parameters and possible offset
Ab = zeros([6,         NChaPTRec+calibrationOffset],'single');%6 motion parameters and 1 offset

if ~useRealImage
    Af = Af + 1i*Af;
    Ab = Ab + 1i*Ab;
end