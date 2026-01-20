
function [parPT] = assignPTPreCalibration(parPT, calibName, forInit)

if nargin<3 || isempty(forInit);forInit=0;end

parPT.Calibration.externalFit=[];
parPT.Calibration.externalFit.calibName=calibName;
if existsFileVar(parPT.Calibration.externalFit.calibName,'PT')==2
    ss = load(parPT.Calibration.externalFit.calibName);
    parPT.Calibration.externalFit.Ab=ss.PT.Ab;
    parPT.Calibration.externalFit.MT=ss.PT.Geom.MT;
    parPT.Calibration.externalFit.NX=ss.PT.Geom.NX;
    if isfield(ss.PT,'preProcessing');parPT.Calibration.externalFit.preProcessing=ss.PT.preProcessing;end
    parPT.Calibration.externalFit.preProcessing.preProcessName = calibName;
end

parPT.Calibration.forInit = forInit;