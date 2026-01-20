
function [rec] = underSampleRec(rec, R, type, dimUS, reOrderSamples, shotsToKeep)

if nargin<2 || isempty(R); R=1;end
if nargin<3 || isempty(type); type = 'regular';end%'incoherent'
if nargin<4 || isempty(dimUS); dimUS = 2;end
if nargin<5 || isempty(reOrderSamples); reOrderSamples = 1;end
if nargin<6 || isempty(shotsToKeep); shotsToKeep = [];end%specic control over NShots

%% SET PARAMETERS
TFE = rec.Par.Labels.TFEfactor;
isShotBased = TFE < length(rec.Assign.z{2});
NShots = length(rec.Assign.z{2}) / TFE;

NY=size(rec.y);

if strcmp(type,'regular')
   if isShotBased
       if NY(2)==TFE;dimUS=3;dimFS=2;else;dimUS=2;dimFS=3;end
   end
   if ~isempty(shotsToKeep)
       linesToKeep = shotsToKeep;
   else
       linesToKeep = round(NY(dimUS)/R);
   end
   fracToKeep = linesToKeep/NY(dimUS);
else 
   if ~isempty(shotsToKeep)
       fracToKeep = shotsToKeep/NShots;
   else
       fracToKeep = 1/R;
   end
end

%% EXTRACT SAMPLES
%%% Compute indices
NSampOld = length(rec.Assign.z{2});
if strcmp(type,'incoherent')
    endSamp = round(fracToKeep*NSampOld);
    if isShotBased;  endSamp = floor(endSamp/rec.Par.Labels.TFEfactor)*TFE;end
    id = 1:endSamp;
    
elseif strcmp(type,'regular')
    %id = mod(rec.Assign.z{dim},R)==0;
    kIndex = PESamplesFromRecPTTest(rec);
    id = mod(kIndex(:,dimUS-1),R)==0;
else
   error('removeKSpaceSamples:: Type not implemented.'); 
end

%%% Extract
for i=2:3; rec.Assign.z{i} = dynInd(rec.Assign.z{i},id,2);end

%% ASSIGNING THE DATA
yOld = rec.y;

%%% Move to k-space
for n=1:3
    yOld=fftshiftGPU(yOld,n);
    yOld=fftGPU(yOld,n)/NY(n);
    yOld=fftshiftGPU(yOld,n);%fftshift since in solveXT there is an iffthift on timeIndex
end
yOld = resSub(yOld, 2:3);

%%% Data to set to 0
[~, idx, hitMat, ~] = PESamplesFromRecPTTest(rec);
rec.y = [];

yNew = zerosL(yOld);
yNew = dynInd(yNew, idx, 2,dynInd(yOld,idx,2));%Zero-filled
yNew = resSub(yNew,2,NY(2:3));

%%% Back to image domain
for n=1:3
    yNew=ifftshiftGPU(yNew,n);
    yNew=ifftGPU(yNew,n)*NY(n);
    yNew=ifftshiftGPU(yNew,n);%fftshift since in solveXT there is an iffthift on timeIndex
end

%%% Store
rec.y=yNew;

%% RE-ORDER
if strcmp(type,'regular') && reOrderSamples
    fprintf('Re-ordering sampling for regular under-sampling.\n')
    zNew = rec.Assign.z;
    [~,idx1] = sort(rec.Assign.z{dimUS},'ascend'); 
    for ii=2:3;zNew{ii} = zNew{ii}(idx1);end
    for i=unique(sort(zNew{dimUS},'ascend'))
        idTemp = zNew{dimUS}==i;
        zTemp = zNew;
        for ii=2:3;zTemp{ii} = zTemp{ii}(idTemp);end
        
        [~,idTemp2] = sort(zTemp{dimFS},'ascend');
        for ii=2:3;zTemp{ii} = zTemp{ii}(idTemp2);end
        
        for ii=2:3;zNew{ii}(idTemp) = zTemp{ii};end
    end
    rec.Assign.z = zNew;
    [~, idx, hitMat, timeMat] = PESamplesFromRecPTTest(rec);
    figure; imshow(timeMat,[]);title('removeKSpaceSamples:: Re-ordering samples to get linear sampling');
    
    %Take only subset of requested
    if linesToKeep * size(rec.y,dimFS)<length(rec.Assign.z{2})
        for ii=2:3;rec.Assign.z{ii} = rec.Assign.z{ii}(1:(linesToKeep * size(rec.y,dimFS)));end
        [~, idx, hitMat, timeMat] = PESamplesFromRecPTTest(rec);
    end
end
NSamNew = length(rec.Assign.z{2});

%%% Change hyper-parameters
if ~isShotBased
    rec.Par.Labels.TFEfactor = length(rec.Assign.z{2}); 
else
    rec.Par.Labels.NShots = length(rec.Assign.z{2})/TFE;
end

%% REPORT
%%% Print
fprintf('removeKSpaceSamples:: %.0f%% kept of samples: from %d to %d samples.\n', 100*fracToKeep, NSampOld, NSamNew);
if isShotBased;fprintf('removeKSpaceSamples:: %d/%d shots kept\n', NSamNew/TFE, NSampOld/TFE);end

%%% Plot
plotND([],RSOS(yNew), defRange(RSOS(yNew)),[],0,[],rec.Par.Mine.APhiRec,[],[],[],100,sprintf('R=%d, %s undersampling',R,type));
figure(101);
imshow(hitMat,[]);title('New samplig mask in PE plane')

