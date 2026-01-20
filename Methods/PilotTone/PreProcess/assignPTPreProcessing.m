
function [parPT] = assignPTPreProcessing(parPT, preProcessName)

if existsFileVar(preProcessName,'PT')==2
    ss = load(preProcessName);
    if isfield(ss.PT,'preProcessing');parPT.signalUsage.preProcessing=ss.PT.preProcessing;end
    parPT.signalUsage.preProcessing.preProcessName = preProcessName;
end