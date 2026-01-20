

function [WeFinal, outlWeFinal] = getPTWeighting_v3(pTime, mSt, NStates, percRobustShot, enerRobustShot, useSVD, nComp, deb, removeNaNs)

%GETPTWEIGHTING finds outlier shots based on the acquired Pilot Tone (PT) signal
%   [WE, OUTLWE] = GETPTWEIGHTING(PTIME, MST, {NSTATES}, {PERCROBUSTSHOT},{ENERGYROBUSTSHOT}, {USESVD}, {NCOMP})
%   * PTIME the pilot tone signal over time with the channels in the rows and time along the columns.
%   * MST the indices defining the grouping.
%   * {NSTATES} the motion states. Note this is only used when there is an elliptical shutter.
%   * {PERCROBUSTSHOT} Percentiles for robust computation of expected inter-shot dispersion.
%   * {ENERGYROBUSTSHOT} Ratio of the error for acceptance of shots.
%   * {USESVD} whether to use the Singular vectors to detect the outliers on.
%   * {NCOMP} the number of Singular vectors to remain for the analysis.
%   * {DEB} is a debug flag.
%   ** WE the soft weights.
%   ** OUTLWE indicating which state is a detected outlier.
%
%   Yannick Brackenier 2023-07-04

if nargin < 3 || isempty(NStates); NStates=length(unique(mSt));end %Number of states
if nargin < 4 || isempty(percRobustShot); percRobustShot=[0.125 0.25];end %Percentile for estimating sigma based on MAD estimator
if nargin < 5 || isempty(enerRobustShot); enerRobustShot=.95;end %ZScore above which to categorise shots as outloers based on their standard deviation
if nargin < 6 || isempty(useSVD); useSVD=1;end %Whether to use the singular components to do the outlier detection or the raw coil signals.
if nargin < 7 || isempty(nComp); nComp=1;end %Number of singular components to keep
if nargin < 8 || isempty(deb); deb=0;end %Number of singular components to keep
if nargin < 9 || isempty(removeNaNs);removeNaNs = 0;end%Don't remove NaNs (0) / Fill with zeros (1) / Remove elements from pGrouped (2)

%%% ENSURE CONSISTENT MST
%[mSt, NStates] = sortmSt(mSt);

%%% PARAMETERS
N = size(pTime);
NCha = N(1); %NSamples = N(2);
nComp = min(nComp,NCha);
percUse=(0.15:0.05:0.85)*100;

%%% CHANNEL COMPRESSION
pTimeProc = pTime;
if deb>=2; visPTSignal (pTimeProc, [], [], [], 1, [], mSt, 100, 'Original signal',[],[],[],[],1);end

if useSVD==1 %Use singular vectors
    [uu, ss, v] = svd(pTimeProc,"econ");
    pTimeProc = v(:,1:nComp).';v=[];
elseif useSVD==2
    pTimeProc = pTimeProc(1:nComp,:);
    pTimeProc = pTimeProc./sqrt(normm(pTimeProc,[],2));
else %Use channels with biggest amplitudes
    v = multDimMea( abs(pTimeProc),2);
    [~,idxMax] = sort(v,'descend');v=[];
    pTimeProc = dynInd(pTimeProc,idxMax(1:nComp),1);
end

%%% GLOBAL DE-MEAN
pTimeProc = pTimeProc - multDimMea(pTimeProc,2);
if deb>=2; visPTSignal(pTimeProc, [], [], [], 1, [], mSt, 101, 'Compressed signal',[],[],[],[],1);end

%%% NORMALISE TO FIRST COMPONENT
if nComp>1
    pTimeProc(2:end,:) = pTimeProc(2:end,:)./sqrt(normm(pTimeProc(2:end,:),[],2)) .* sqrt(normm(pTimeProc(1,:),[],2));
    if deb>=2; visPTSignal (pTimeProc, [], [], [], 1, [], mSt, 102, 'Compressed + normalised signal',[],[],[],[],1);end
end

%%% SHOT HANDLING
pTimeProcGr = groupPT(pTimeProc, mSt(mSt<=NStates), NStates);
idxFilled = ~isnan(pTimeProcGr(1,:));
for s=1:NStates; pTimeProc(:,mSt==s) = abs(pTimeProc(:,mSt==s) - pTimeProcGr(:,s) );  end
if deb>=2; visPTSignal (pTimeProc, [], [], [], 1, [], mSt, 103, 'Compressed + normalised + de-meaned signal',[],[],[],[],1);end

%%% AVERAGE
pTimeProc = multDimMea(pTimeProc,1);
if deb>=2; visPTSignal(pTimeProc, [], [], [], 1, [], mSt, 104, 'Compressed + normalised + de-meaned + averaged signal',[],[],[],[],1);end

%%% WEIGHT CALCULATION
WeO=zeros([NStates length(percUse)],'like',pTimeProc);
for s=1:NStates;WeO(s,:)=prctile(pTimeProc(mSt==s),percUse);end
    WeO = dynInd(WeO,idxFilled,1);
WeO=permute(WeO,[1 3 2]);           
Westd=diff(prctile(WeO,percRobustShot*100,1),1,1)./diff(norminv(percRobustShot)); 
Wmea=prctile(WeO,mean(percRobustShot)*100,1);
Wmea=Wmea+Westd*sqrt(2)*erfcinv(2*mean(percRobustShot));%sqrt(2)*erfcinv(2*mean(parXT.percRobustShot))=norminv(mean(parXT.percRobustShot))
We=bsxfun(@minus, WeO,Wmea);
We=bsxfun(@rdivide, We, Westd);
We=mean(We,3);         
We=min((1-normcdf(We))/((1-enerRobustShot)/multDimSum(idxFilled)),1);  

%%% OUTLIER DETECTION
outlWe=We<1; 

%%% FILL
WeFinal = NaN([1 NStates],'like',We);
WeFinal(idxFilled) = We;
outlWeFinal = single(NaN([1 NStates]));
if isa(We, 'gpuArray'); outlWeFinal=gpuArray(outlWeFinal);end
outlWeFinal(idxFilled) = single(outlWe);

%%% REMOVE EMPTY STATES
if removeNaNs==1%Fill with zero
    WeFinal = dynInd(WeFinal,~idxFilled,2,1);
    outlWeFinal = dynInd(outlWeFinal,~idxFilled,2,0);
elseif removeNaNs==2%Remove NaNs
    WeFinal = dynInd(WeFinal,idxFilled,2);
    outlWeFinal = dynInd(outlWeFinal,idxFilled,2);
end

%%% PLOT
if deb>=1
    visPTSignal (pTime, [], [], [], 1, outlWeFinal, mSt, 105, 'Outlier detection',[],[],[],[],1);
    visResiduals(We,outlWe,gather(1:multDimSum(idxFilled)),[],[],0, [],[],106,'Weights');
end




