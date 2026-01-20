
function [WeFinal, outlWeFinal] = computeWe(yX,xRes,E,EH,iSt,l,n,NYres,NSweeps,parXT,mStSweeps,leff, Re)

if nargin<13 || isempty(Re);Re=[];end

%COMPUTE RESIDUALS IF NOT PROVIDED
if isempty(Re)
    [Re{l}(:,n),~,fullResid]=computeEnergy(yX,xRes,E,[],EH,[],[],[],[],1);Re{l}(:,n)=Re{l}(:,n)/numel(yX);
end
Residuals=Re{l}(:,n);Residuals(iSt{l})=Residuals;
Residuals=reshape(Residuals,NYres([1:2 5]));

mStSort=mStSweeps;mStSort(iSt{l})=mStSort;%YB: now in mStSort the samples correspond to the right sweep
Sweeps=reshape(mStSort,NYres([1:2 5]));%YB: So if we reshape we have the PE plane with the sweep number per TR         

%OUTLIER DETECTION
if leff>0;percUse=(0.15:0.05:0.85)*100;
else; percUse=(0.3:0.05:0.7)*100;
end
WeO=zeros([NSweeps length(percUse)],'like',Residuals);
for s=1:NSweeps;WeO(s,:)=prctile(log(Residuals(Sweeps==s)),percUse);end
    idxFilledStates = ~isnan(WeO(:,1,1));
    WeO = dynInd(WeO,idxFilledStates,1);
    NSweepsTemp = multDimSum(idxFilledStates);
WeO=permute(WeO,[1 3 2]);           
Westd=diff(prctile(WeO,parXT.percRobustShot*100,1),1,1)./diff(norminv(parXT.percRobustShot)); 
Wmea=prctile(WeO,mean(parXT.percRobustShot)*100,1);
Wmea=Wmea+Westd*sqrt(2)*erfcinv(2*mean(parXT.percRobustShot));%sqrt(2)*erfcinv(2*mean(parXT.percRobustShot))=norminv(mean(parXT.percRobustShot))
We=bsxfun(@minus, WeO,Wmea);%YB changed for version R2015
We=bsxfun(@rdivide, We, Westd);%YB changed for version R2015
We=mean(We,3);         
We=min((1-normcdf(We))/((1-parXT.enerRobustShot)/NSweepsTemp),1);%PTHandling: even if WeO==NaN, that is fine as this will be set to 1 in this line        

%Fill the states that contained samples
WeFinal = ones([NSweeps 1],'like', We);
WeFinal = dynInd(WeFinal, idxFilledStates,1, We);
outlWeFinal = WeFinal<1;

