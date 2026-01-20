

function [TResid, limT, TPT] = clusteringDifferenceRes(PT, T, stateSample, MS)

%%% COMPUTE MOTION TRACE FROM PILOT TONE
TPT = pilotTonePrediction(PT.Ab, PT.pTimeRes, 'PT2T','backward', PT.NChaPTRecCurr, 6);

%%% REMOVE ALIGNED-SENSE MOTION ESTIMATES (AT TR LEVEL)
TResid = TPT;
for state = 1:max(stateSample)
    stateIdx = stateSample==state;
    stateError = dynInd(TResid, stateIdx, 5) - dynInd(T,state,5);
    TResid = dynInd(TResid, stateIdx, 5,stateError );
end

%%% COMPUTE MOTION RANGE
limT=getTLims(TResid,MS);
