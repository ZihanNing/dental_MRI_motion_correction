
function [outlPTNew] = plotPTSignalsToCompare_v2(pGT, p1, p2, We,numFig, voxSize)
%just a temporary function to see how motion/PT signals differ when using
%the forward and backward model

if nargin < 5 || isempty(numFig);numFig=90;end
if nargin < 6 || isempty(voxSize);voxSize=[1 1 1];end

%%% REHAPE IF p are actual motion parameters
N = size(pGT);
isMotion = length(N)>2;

if ~isempty(We); idxGoodShot = We==1;else;idxGoodShot = ones([1 N(5)]);We = onesL(idxGoodShot);end

if isMotion
    %Permute so NCha x NStates
    pGT=permute(pGT, [6 5 1:4]); 
    p1=permute(p1, [6 5 1:4]); 
    p2=permute(p2, [6 5 1:4]); 
    %Make rotation angles in degrees
    pGT(4:6,:) = 180/pi*pGT(4:6,:);
    p1(4:6,:) = 180/pi*p1(4:6,:);
    p2(4:6,:) = 180/pi*p2(4:6,:);
    
    pGT(1:3,:) = voxSize(:).*pGT(1:3,:);
    p1(1:3,:) = voxSize(:).*p1(1:3,:);
    p2(1:3,:) = voxSize(:).*p2(1:3,:);
    
else
    pGT=abs(pGT);
    p1=abs(p1);
    p2=abs(p2);
end
NSamples= size(pGT,2);
NCha= size(pGT,1);

%%%FIGURE
h  = createFig(numFig);
i=0;
idxTemp = idxGoodShot==i;

%%% first set 
subplot(3,2,1); plot(1:multDimSum(idxTemp), dynInd(pGT,idxTemp,2))
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('OUTLIERS: alignedSENSE parameters')
axis tight

subplot(3,2,3); plot(1:multDimSum(idxTemp), dynInd(p2,idxTemp,2))
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('Prediction - forward model')
axis tight

subplot(3,2,5); plot(1:multDimSum(idxTemp), abs(dynInd(pGT,idxTemp,2)-dynInd(p2,idxTemp,2)) )
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('Difference with measurement - forward model')
axis tight

%%% second set
i=1;
idxTemp = idxGoodShot==i;
subplot(3,2,2); plot(1:multDimSum(idxTemp), (dynInd(pGT,idxTemp,2)))
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('RELIABLE SHOTS: alignedSENSE parameters')
axis tight

subplot(3,2,4); plot(1:multDimSum(idxTemp), (dynInd(p2,idxTemp,2)))
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('Prediction - backward model')
axis tight

subplot(3,2,6); plot(1:multDimSum(idxTemp), abs(dynInd(pGT,idxTemp,2)-dynInd(p2,idxTemp,2) ) )
xlabel('Motion states')
if ~isMotion; ylabel('PT magnitude');else; ylabel('Motion parameters');end
title('Difference with measurement - backward model')
axis tight

set(h,'color','w','Position',get(0,'ScreenSize'));
if isMotion; titleName='Motion';else;  titleName='PT';end
sgtitle(sprintf('Comparing  %s parameters after calibration',titleName));

if nargout>0
    res = abs(pGT(1:3,:)-p2(1:3,:));
    outlPTNewTran = multDimMea(res,1)>multDimMea(voxSize)/10;

    res = abs(pGT(4:6,:)-p2(4:6,:));
    radiusHead = 95;%mm
    outlPTNewRot = multDimMea(res,1)>asind(multDimMea(voxSize)/10/radiusHead);

    outlPTOrig = We<1;
    outlPTNew = outlPTOrig | outlPTNewTran | outlPTNewRot;
end
