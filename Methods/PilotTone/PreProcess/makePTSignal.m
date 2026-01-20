

function pTimeRes = makePTSignal(pSliceRes,NYres,timeIndexRes,useRealImag,NChaPTRec)

pTimeRes = permute(pSliceRes,[1 2 5 3 4]);%5 is for the multiple repeats/shots
pTimeRes = reshape(pTimeRes,[prod(NYres([1:2 5])) 1 1 size(pTimeRes,5)]);%pTimeRes = resPop(pTimeRes ,1:3,[],1);

[timeIndexResSorted,idxx]=sort(timeIndexRes(:),'ascend');idxx(timeIndexResSorted==0)=[];
pTimeRes = dynInd( pTimeRes,idxx,1);%Within each grouping, samples are not sorted in temporal order, so plot might look weird. Does not matter since you average samples within the group
pTimeRes = permute(pTimeRes , [4 1 2:3]);%Make PT signal as NCha x NStates           
pTimeRes = extractPTType(pTimeRes,useRealImag, NChaPTRec);
        