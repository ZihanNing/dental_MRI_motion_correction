

function [pTimeRes, stateSampleRes, stateSampleIdx] = getPTSignalRes(PT,parXT,timeIndexRes,NYres,stateSample)

%GETPTSIGNALRES extracts the PT signal at a given resolution.
%   [PTIMERES,STATESAMPLERES,STATESAMPLEIDX]=GETPTSIGNALRES(PT,PARXT,TIMEINDEXRES,NYRES,STATESAMPLE)
%   * PT is PT structure in rec that contains all the PT data. 
%   * PARPT is PT structure in rec.Alg that contains all the PT processing information. 
%   * TIMEINDEXRES is 2D array containing the time index of the samples at the given resolution. 
%   * NYRES is the k-space data at the resolution. 
%   * STATESAMPLE is the array of size [1 NSamples] that contains the state of each sample at the native resolution. 
%   ** PTIMERES is sorted PT signal at that resolution.
%   ** STATESAMPLERES is stateSample array of the corresponding resolution.
%   ** STATESAMPLEIDX is an array with indices of the array stateSamples whose elements are present at the current resolution.
%   ** PTIMERES is sorted PT signal at that resolution.
%
%   Yannick Brackenier 2023-03-17

PT.pTimeRes = permute(PT.pSliceRes,[1 2 5 3 4]);%5 is for the multiple repeats/shots
PT.pTimeRes = reshape(PT.pTimeRes,[prod(NYres([1:2 5])) 1 1 size(PT.pTimeRes,5)]);%PT.pTimeRes = resPop(PT.pTimeRes ,1:3,[],1);

[timeIndexResSorted,idxx]=sort(timeIndexRes(:));idxx(timeIndexResSorted==0)=[];
PT.pTimeRes = dynInd(PT.pTimeRes,idxx,1);%Within each grouping, samples are not sorted in temporal order, so plot might look weird. Does not matter since you average samples within the group
PT.pTimeRes = permute(PT.pTimeRes , [4 1 2:3]);%Make PT signal as NCha x NStates           
PT.pTimeRes = extractPTType(PT.pTimeRes, PT.useRealImag, PT.NChaPTRecCurr);
      
pTimeRes = PT.pTimeRes;
timeIndexRes = timeIndexRes(timeIndexRes~=0);
stateSampleIdx = sort(timeIndexRes);
stateSampleRes = stateSample(stateSampleIdx);        
        












%         
%     E.mSt=timeIndexRes;  
%     E.mSt(timeIndexRes~=0)=stateSample(timeIndexRes(timeIndexRes~=0)); %YB: Now mSt a 2D grid with the same states having the value of the state number. Also note that timeIndexRes contains indices of the samples, so biggere than numel(timeIndexRes)
%     
%     if estT(l) || ~isfield(EH,'We');EH.We=[];end    
%     if parXT.fractionOrder~=0;GRes=buildFilter([NYres(1:2) length(E.nF)],'FractionalFiniteDiscreteIso',NX./NXres,gpuIn,parXT.fractionOrder);else; GRes=[];end%For fractional finite difference motion estimation    %YB: was size NXres but changed to NYres since it needs to cover k-space (matters for accelerated scans)       
%     fprintf('Building harmonic sampling structures.\n');
%     [ySt,FSt,GRes,E.iSt,E.nSt,E.mSt,nSa]=buildHarmonicSampling(yRes,GRes,E.mSt,NStates,NXres(1:2));iSt{l}=cat(1,E.iSt{:});yRes=gather(yRes);%Samples as cell arrays
%     if any(nSa(1:NStates)==0)
%         %nowithin=1;fprintf('Probably not DISORDER, returning without performing corrections\n');return;   
%         warning('Some shots have no samples.');
%     end
   
