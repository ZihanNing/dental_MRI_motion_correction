
function [recS] = rec2recS(rec)

recS = rec;

recS.y=rec.ACS;
recS = rmfield(recS,'ACS');

recS.Par.Mine.APhiRec = recS.Par.Mine.APhiACS;

recS.Enc.AcqVoxelSize = MT2MS(recS.Par.Mine.APhiRec);

recS.Enc.FOVSize = multDimSize(recS.y,1:3); 
recS.Enc.AcqSize =  multDimSize(recS.y,1:3); 
recS.Enc.AcqSize(1) = 2*recS.Enc.AcqSize(1);

recS.Enc.UnderSampling = [];
recS.Enc.UnderSampling.Flag = 0;

recS.Assign = [];

recS.Enc.DISORDER.Flag = 0;
recS.Enc.kRange = [];
