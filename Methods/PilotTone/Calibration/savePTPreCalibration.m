
function [] = savePTPreCalibration(calibName, Ab, MT, NX, preProcessing )

if nargin<5 || isempty(preProcessing);preProcessing=[];end

PT=[];

PT.Ab=Ab;
PT.Geom.MT=MT;
PT.Geom.NX=NX;
if ~isempty(preProcessing);PT.preProcessing=preProcessing;end

%save
save(calibName,'PT');