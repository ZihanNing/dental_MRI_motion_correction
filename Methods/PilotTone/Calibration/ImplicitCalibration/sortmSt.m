

function [mStSorted, NStates] = sortmSt(mSt,NStates)

%SORTMST sorts the array with the motion state indices of the samples.
%   [MSTSORTED]=SORTMST(MST,{NSTATES})
%   * MST is array of size [1 NSamples] containing the motion state index of each sample. 
%   * {NSTATES} is the total number of states. Defaults to the maximum values in mSt. Larger values can exist (e.g. in case of sequential sampling).
%   ** MSTSORTED is sorted version with increasing indices.
%
%   Yannick Brackenier 2023-03-17

uniqueSates = sort(unique(mSt));
newStates = 1:length(uniqueSates);
if nargin<2 || isempty(NStates);NStates=length(uniqueSates);end
assert(NStates==length(uniqueSates),'sortmSt:: NStates specified does not agree with the provided mSt.')

mStSorted=zerosL(mSt);
for i=1:NStates
    idx = mSt==uniqueSates(i);
    mStSorted(idx)=newStates(i);
end
assert(~any(mStSorted==0),'sortmSt:: State 0 not allowed.')
