
function []= visMotionFromMat(fileName)

load(strcat(fileName,'.mat'));

voxSize = MT2MS(MotionInfo.Par.Mine.APhiRec);
voxSize = voxSize(MotionInfo.Par.Mine.permuteHist{end}(1:3));

%DISORDER motion trace
MotionInfo.T = T;
time = multDimMea(MotionInfo.timeState,2:ndims(MotionInfo.timeState));
limTraces=[];
visMotion(MotionInfo,voxSize,time,2,[],[],limTraces,[],1);
