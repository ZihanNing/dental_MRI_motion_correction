
function pGrouped = groupPT(p, stateSample, NStates, removeNaNs)

%GROUPPT groups Pilot Tone (PT) signal into groups based on a state (group) index.
%   [PGROUPED]=GROUPPT(P, STATESAMPLE, {NSTATES})
%   * P the PT signal as an array with size [NCha x NSamples].
%   * STATESAMPLE the state (group) index for each sample. Array with size [1 x NSamples].
%   * {NSTATES} the number of states (groups). Defaults to max(stateSample).
%   * {REMOVENANS} a flag to remove or replace the NaNs associated with empty states.
%   ** PGROUPED is the grouped PT signal as an array with size [NCha x NStates].
%
%   Yannick Brackenier 2023-04-05

if isa(stateSample, 'gpuArray'); stateSample=gather(stateSample);end

NStateCurrent = max(stateSample(:));
if nargin<3 || isempty(NStates);NStates = NStateCurrent;end
if nargin<4 || isempty(removeNaNs);removeNaNs = 0;end%Don't remove NaNs (0) / Fill with zeros (1) / Remove elements from pGrouped (2)

N = size(p);
NCha = N(1); %NSamples = N(2);

%%% CREATE ARRAY
pGrouped = NaN([NCha, NStateCurrent],'like', p);%Initial size [NCha, NStateCurrent] since NStateCurrent will come out of regionprops.
isComplex = ~isreal(p);

%%% PERFORM GROUPING
for cha=1:NCha %Each coil individually
    %Real component
    ttReal = regionprops(stateSample,dynInd(real(p),cha,1),'MeanIntensity');%Can probably made faster with accumarray.m
    ttReal = {ttReal.MeanIntensity};
    ttReal = cat(2,ttReal{:});
    
    %Imaginary component
    if ~isComplex
        pGrouped = dynInd(pGrouped,cha,1,ttReal); 
    else
        ttImag = regionprops(stateSample,dynInd(imag(p),cha,1),'MeanIntensity');
        ttImag = {ttImag.MeanIntensity};
        ttImag = cat(2,ttImag{:});
        pGrouped = dynInd(pGrouped,cha,1,ttReal + 1i*ttImag); 
    end
end

%%% EXTEND ARRAY FROM NSTATESCURRENT TO NSTATES
if isComplex; valToFill = NaN + 1i*NaN; else; valToFill=NaN;end
pGrouped = dynInd(pGrouped,NStateCurrent+1:NStates,2,valToFill);%Extend with NaNs at the end where you might have zero-samples shots

%%% REMOVE EMPTY STATES
idxEmpty = isnan(pGrouped(1,:));
if removeNaNs==1%Fill with zero
    pGrouped = dynInd(pGrouped,idxEmpty,2,0);
elseif removeNaNs==2%Remove NaNs
    pGrouped = dynInd(pGrouped,~idxEmpty,2);
end


