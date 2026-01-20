
function [yFit, x] = curveFit(y, forwardModel, initGuess, paramBounds, eqCon, nonEqCon, nonLocalCon, dimFit, sepFit)

if nargin <3 || isempty(initGuess);initGuess=[];end
if nargin <4 || isempty(paramBounds);paramBounds=[];end
if nargin <5 || isempty(eqCon);eqCon=[];end
if nargin <6 || isempty(nonEqCon);nonEqCon=[];end
if nargin <7 || isempty(nonLocalCon);nonLocalCon=[];end
if nargin <8 || isempty(dimFit);dimFit=numDims(y);end
if nargin <9 || isempty(sepFit);sepFit=0;end

%%% RESHAPE
% perm = 1:numDims(y);perm(dimFit)=[];perm=[perm dimFit];
% y = permute(y,perm);
% NYPerm = size(y);
% y = resSub(y,1:(length(NYPerm)-1));
% NG = size(y,2);%Number of elements along dimension on which to fit curve
% NS = size(y,1);%Number of instances (separable fits)

%%% SET OPTIMISATION PARAMETERS
if isstruct(nonEqCon);A = nonEqCon.A; b = nonEqCon.b;else; A=[];b=[];end
if isstruct(eqCon);A = eqCon.A; b = eqCon.b;else; Aeq=[];beq=[];end

if size(paramBounds, 1)>0; lowerBound = dynInd(paramBounds,1,1);else; lowerBound = [];end
if size(paramBounds, 1)>1; upperBound = dynInd(paramBounds,2,1);else; upperBound = [];end
if isempty(initGuess) 
    if ~(isempty(lowerBound) || isempty(upperBound))
        initGuess = (lowerBound + upperBound)/2;
    else
        initGuess = ones([ 1 funHandleInputSize(forwardModel)]);
        warning('curveFit:: No initial values set. All set to 1.');
    end 
end
%if size(initGuess,1)~=size(y,1); initGuess = cat(1,initGuess,repmat(initGuess(1,:), [size(y,1)-size(initGuess,1) ,1] ) );end
    
optionsCon = optimoptions('fmincon','Algorithm','sqp','Display','none','Diagnostics','off');
%optionsCon = optimoptions('fmincon','Display','iter','Diagnostics','off');
optionsCon.MaxFunctionEvaluations=5*1e3;
optionsCon.StepTolerance=1e-8;
optionsCon.OptimalityTolerance=1e-8;
%optionsUnCon = optimset('fminsearch','Algorithm','sqp');
% optionsUnCon = [];    
% function f = lossFun(x); f = double(normm(y,forwardModel(x))); end

%%% SET TO DOUBLE FOR FMINCON EXECUTION
y=double(y);
initGuess=double(initGuess);

%%% FIT DATA
if sepFit
    x = zerosL(initGuess);
    blSz = 1;
    for i=1:blSz:size(x,1)
        %i
        idx=i:min(i+blSz-1,size(x,1));
        lossFun = @(x) double(normm(y(idx,:),forwardModel(x))); 
        [xTemp,~,exitflag,~] = fmincon(lossFun,initGuess(idx,:),A,b,Aeq,beq,lowerBound,upperBound,nonLocalCon,optionsCon);
        x = dynInd(x,idx,1,xTemp);
        if ~(exitflag>0); fprintf('curveFit:: fmincon not converged.');end
    end
else
    lossFun = @(x) double(normm(y,forwardModel(x))); 
    [x,~,exitflag,~] = fmincon(lossFun,initGuess,A,b,Aeq,beq,double(lowerBound(:)),double(upperBound(:)),nonLocalCon,optionsCon);
    %[x,~,exitflag,~] = fminsearch(lossFun,initGuess,[]);
    if ~(exitflag>0); fprintf('curveFit:: fmincon not converged.');end
    lossFun(initGuess)
    lossFun(x)
    
end
%%% RE-EVALUATE FORWARD MODEL
yFit = forwardModel(x);

%%% SET TO SINGLE 
yFit=single(yFit);
x=single(x);

%%% RESHAPE INTO ORIGINAL SIZE
% yFit = resSub(yFit, 1, NYPerm(1:(length(NYPerm)-1)));
% yFit = ipermute(yFit,perm);

end